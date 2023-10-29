#include <iostream>

#include "a9.h"
#include "Image.h"
#include "basicImageManipulation.h"
#include "matrix.h"
#include "filtering.h"
#include "morphing.h"

using namespace std;

// Write your implementations here, or extend the Makefile if you add source
// files
void someFunction() {
    cout << "ok, that's a function" << endl;
}

float dotIm(const Image &im1, const Image &im2, const bool &zeroClamp){
    float val =0.0;
    for (int x=0; x<im1.width(); x++){
        for (int y=0; y<im1.height(); y++){
            for (int z=0; z<im1.channels(); z++){
                val = val + im1(x,y,z)*im2(x,y,z);
            }
        }
    }
    return val;
}
Image applyKernel(const Image &im, Filter &filt){
    return filt.convolve(im);
}
Image applyConjugatedKernel(const Image &im, Filter &filt){
    return filt.transpose().convolve(im);
}
Image computeResidual(const Image &x, const Image &y, Filter &kernel){
    return y - kernel.convolve(x);
}


Filter laplacian(){
    return Filter(vector<float>({0,-1,0,-1,4,-1,0,-1,0}),3,3);
}
Image applyLapacian(const Image &im){
    return laplacian().convolve(im);
}
Image applyRegularizedOperator(const Image &x, Filter &kernel, const float &lamb){
    return kernel.transpose().convolve(kernel.convolve(x)) + lamb*laplacian().convolve(x);
}

Image deconvCG_reg(const Image &y, Filter &kernel, const float &lamb, const int &niter){
    Image x = Image(y.width(),y.height(),y.channels());
    Image r = applyConjugatedKernel(y,kernel) - applyRegularizedOperator(x,kernel,lamb);
    Image d = r;
    Image Ad =x;
    float r_dot, alpha, beta ;
    for (int i=0; i<niter; i++){
            Ad = applyRegularizedOperator(d,kernel,lamb);
            r_dot = dotIm(r,r);
            alpha = r_dot/dotIm(d,Ad);  
            x = x + alpha*d;
            r = r - alpha*Ad;
            beta = dotIm(r,r)/r_dot;
            d = r + beta*d;
        }
    return x;
}
Image poisson(const Image &bg, const Image &fg, const Image &mask, const int niter){
    Image x = (1-mask) * bg;
    Filter lap = laplacian();
    Image b = lap.convolve(fg);
    Image r = mask*b;
    Image Ar =Image(bg.width(),bg.height(),bg.channels()) ;
    float alpha ;
    for (int i=0; i<niter; i++){
        Ar = lap.convolve(r);
        // cout <<r.width()<<' '<<r.height()<<' '<<r.channels()<<endl;
        alpha = dotIm(r,r)/dotIm(r,Ar);  
        x = x + alpha*r;
        r = mask*(b - lap.convolve(x));
    }
    return x;
}

Image poissonCG(const Image &bg, const Image &fg, const Image &mask, const int niter){
    Image x = (1-mask) * bg;
    Filter lap = laplacian();
    Image r = mask*(lap.convolve(fg) - lap.convolve(x));
    Image d = r;
    Image Ad =x;
    float r_dot, alpha, beta ;
    for (int i=0; i<niter; i++){
        Ad = lap.convolve(d);
        // cout <<r.width()<<' '<<r.height()<<' '<<r.channels()<<endl;
        r_dot = dotIm(r,r);
        alpha = r_dot/dotIm(d,Ad);  
        x = x + alpha*d;
        r = mask*(r - alpha*Ad);
        beta = dotIm(r,r)/r_dot;
        d = r + beta*d;
    }
    return x;
}


vector<Image> channel_split(const Image &im){
    Image imr = Image(im.width(),im.height(),1);
    Image img = Image(im.width(),im.height(),1);
    Image imb = Image(im.width(),im.height(),1);
    vector<Image> res = vector<Image>({imr,img,imb});
    for (int z=0;z<im.channels();z++){
        for (int x=0;x<im.width();x++){
            for (int y=0; y<im.height();y++){
                res[z](x,y) = im(x,y,z);
            }
        }
    }
}

Image channel_combine(const vector<Image> &IM){
    Image res = Image(IM[0].width(),IM[1].height(),3);
    for (int z=0;z<3;z++){
        for (int x=0;x<res.width();x++){
            for (int y=0; y<res.height();y++){
                res(x,y,z) = IM[z](x,y);
            }
        }
    }
}

Image channel_copy(const Image &im){
    Image res = Image(im.width(),im.height(),3);
    for (int x=0; x<im.width(); x++){
        for (int y=0; y<im.height(); y++){
            for (int z=0; z<3; z++){
                res(x,y,z)= im(x,y);
            }
        }
    }
    return res;
}

Image computeTensor(const Image &im, float sigmaG, float factorSigma) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Compute xx/xy/yy Tensor of an image. (stored in that order)
  Image imb = gaussianBlur_separable(im,sigmaG);
  Image imx = gradientX(imb,true);
  Image imy = gradientY(imb,true);
  Image res = Image(im.extent(0),im.extent(1),3);
  for (int x =0; x<res.width();x++){
    for (int y=0; y<res.height();y++){
      res(x,y,0) = powf(imx(x,y),2);
      res(x,y,1) = imx(x,y)*imy(x,y);
      res(x,y,2) = imy(x,y)*imy(x,y);
    }
  }
  Image result = gaussianBlur_separable(res,sigmaG*factorSigma);

  return result;
}

Image inpaint_poisson(const Image &im, const Image &mask){
    Image flat = Image(im.width(), im.height(), im.channels());
    // cout <<sTensor.width()<<sTensor.height()<<flat.channels()<<' '<<mask.width()<<mask.height()<<mask.channels()<<endl;
    return poissonCG(im, flat,mask,200);
}

float solve(int i1,int j1,int i2,int j2, Image &f, Image &T,int KNOWN){
    float sol = 1.0e6;
    if (f(i1,j1)==KNOWN){
        if (f(i2,j2)==KNOWN){
            float r = sqrt(2*(T(i1,j1)*T(i2,j2))*(T(i1,j1)*T(i2,j2)));
            float s = (T(i1,j1)+T(i2,j2)*r)/2;
            if (s>=T(i1,j1) && s>=T(i2,j2)){ 
                sol = s;
            }else{ 
                s += r; if (s>=T(i1,j1) && s>=T(i2,j2)) sol = s; 
            }
        }else{
            sol = 1+T(i1,j1);
        }
    }else if (f(i2,j2)==KNOWN) {
        sol = 1+T(i1,j2);
    }
    return sol;
}

void inpaint_point(int i,int j, const Image &mask, const Image &T, const Image &Tx, const Image &Ty, Image &out, int eps){
    float w, s, Ia, dir, dst, lev;
    Vect2f r, gradI;
    float OUTSIDE =0;
    for (int k=i-eps; k<=i+eps; k++){
        for (int l=j-eps; l<=j+eps; l++){
            if (mask(k,l)==OUTSIDE){
                r = Vect2f(i,j) - Vect2f(k,l); //= vector from (i,j) to (k,l);
                dir = dot(r, Vect2f({Tx(i,j),Ty(i,j)}))/length(r);
                dst = 1/(length(r)*length(r));
                lev = 1/(1+fabs(T(k,l)-T(i,j)));
                w = dir*dst*lev;
                if (mask(k+1,l)!=OUTSIDE && mask(k-1,l)!=OUTSIDE && mask(k,l+1)!=OUTSIDE && mask(k,l-1)!=OUTSIDE){
                    gradI = Vect2f(out(k+1,l)-out(k-1,l),out(k,l+1)-out(k,l-1));
                    Ia += w * (out(k,l) + dot(gradI, r));
                    s += w;
                }
            }
        }
    }
    out(i,j) = Ia/s;
}
float getx(const Image &im, const int x){
    if (x>=im.number_of_elements() || x<0){
        return 1.0f;
    }else{
        return im(x);
    }
}
Image inpaint(const Image &im, const Image &mask){
    vector<Image> img = channel_split(im);
    vector<int> narrow_band; // TODO maybe some other data struct
    int w = im.width();
    int p,i,j;
    Image f = mask; //TODO band image
    Image T = 1.0e06*mask;
    int k[4], l[4];
    int KNOWN =0;
    int INSIDE = 1;
    for (int x=0; x<mask.number_of_elements(); x++){
        if (mask(x)==0 && (getx(mask,x-1)==1 || getx(mask, x+1)==1 || getx(mask,x+w)==1 || getx(mask,x-w)==1)){
            narrow_band.push_back(x);
            f(x) = 0.5;
        }
    }
    while (narrow_band.size()>0){
        p = narrow_band[0]; 
        f(p) = KNOWN;
        i= p%w; 
        j= (p-i)/w;
        k[0] =i-1;k[1]=i;k[2]=i+1;k[3]=i;
        l[0] =j;k[1]=j-1;k[2]=j;k[3]=j+1;
        for (int u=0; u<4;u++){
            if (f(k[u],l[u])!=KNOWN){
                // T(k[u],l[u]) = min({
                //     solve(k[u]-1, l[u], k[u], l[u]-1,f,T),
                //     solve(k[u]+1, l[u], k[u], l[u]-1,f,T),
                //     solve(k[u]-1, l[u], k[u], l[u]+1,f,T),
                //     solve(k[u]+1, l[u], k[u], l[u]+1,f,T)
                // });
                narrow_band.push_back(l[u]*img[0].width()+k[u]); //TODO
            }
        }

    }
    Image Tx = gradientX(T);
    Image Ty = gradientY(T);
    for (int x=0; x<img[0].width(); x++){
        for (int y=0; y<img[0].height(); y++){
            inpaint_point(x,y,mask,T, Tx, Ty, img[0]);
            inpaint_point(x,y,mask,T, Tx, Ty, img[1]);
            inpaint_point(x,y,mask,T, Tx, Ty, img[2]);
        }
    }
    return channel_combine(img);
}