// npr.cpp
// Assignment 8

#include "npr.h"
#include "filtering.h"
#include "matrix.h"
#include <algorithm>
#include <math.h>

using namespace std;

/**************************************************************
 //                       NPR                                //
 *************************************************************/
float randm(){
	return rand()/((float)(RAND_MAX));
}

void brush(Image &im, int x, int y, vector<float> color, const Image &texture) {
	// Draws a brushstroke defined by texture and color at (x,y) in im
	// // --------- HANDOUT  PS12 ------------------------------
	int xs = max(x -texture.width()/2, texture.width()/2), xe=min(x+texture.width()/2, im.width()-texture.width()/2);
	int ys = max(y -texture.height()/2, texture.height()/2), ye=min(y+texture.height()/2, im.height()-texture.height()/2);
	float op;
	int x2,y2;
	for (int x1=0; x1<texture.width(); x1++){
		for (int y1=0; y1<texture.height(); y1++){
			x2 = x + x1 - texture.width()/2;
			y2 = y + y1 - texture.height()/2;
			if (x2>=xs and x2<=xe and y2>=ys and y2<=ye){
				op = texture(x1,y1);
				for (int z=0;z<3;z++){
					im(x2,y2,z) = (1-op)*im(x2,y2,z) + op *color[z];
				}
			}

		}
	}
	return;
}

void singleScalePaint(const Image &im, Image &out, const Image &texture, int size, int N, float noise) {
	// Create painted rendering by splatting brushstrokes at N random locations
	// in your ouptut image
	// // --------- HANDOUT  PS12 ------------------------------
	float factor = min(float(size)/texture.width(),float(size)/texture.height());
	Image tscaled = scaleLin(texture, factor);
	float nfact = 0;
	int x, y;
	vector<float> color = {0,0,0};
	for (int i=0;i<N;i++){
		x = int(randm()*(out.width()-1));
		y = int(randm()*(out.height()-1));
		for (int z=0;z<3;z++){
			nfact = 1 - noise/2 + noise*randm();
			color[z] = nfact*im(x,y,z);
		}
		if (factor<1){
			brush(out,x,y,color,tscaled);
		}else{
			brush(out,x,y,color,texture);
		}
		
	}
	return;
}

void singleScalePaintImportance(const Image &im, const Image &importance,
						Image &out, const Image &texture, int size, int N, float noise) {
	// Create painted rendering but vary the density of the strokes according to
	// an importance map
	// // --------- HANDOUT  PS12 ------------------------------
	float factor = min(float(size)/texture.width(),float(size)/texture.height());
	Image tscaled = texture;
	float nfact = 0;
	int x, y;
	vector<float> color = {0,0,0};
	if (factor<1){
		tscaled = scaleLin(texture, factor);
	}
	float sum=0;
	for (int x3=0;x3<importance.width();x3++){
		for (int y3=0;y3<importance.height();y3++){
			sum = sum + importance(x3,y3);
		}
	}
	int N1 = N/(sum/(importance.width()*importance.height()));
	for (int i=0;i<N1;i++){
		x = randm()*(out.width()-1);
		y = randm()*(out.height()-1);
		for (int z=0;z<3;z++){
			nfact = 1 - noise/2 + noise*randm();
			color[z] = nfact*im(x,y,z);
		}
		if (importance(x,y)>=randm()){
			brush(out,x,y,color,tscaled);
		}
	}

	return;
}

Image sharpnessMap(const Image &im, float sigma) {
	// Calculate sharpness mask 
	// // --------- HANDOUT  PS12 ------------------------------
	vector<Image> lc = lumiChromi(im);
	Image hf = lc[0] - gaussianBlur_separable(lc[0],sigma);
	hf = hf*hf;
	hf = gaussianBlur_separable(hf, 4*sigma);
	hf = hf/hf.max();
	return hf;
}

void painterly(const Image &im, Image &out, const Image &texture, int N, int size, float noise) {
	// Create painterly rendering using a first layer of coarse strokes followed
	// by smaller strokes in high detail areas
	// // --------- HANDOUT  PS12 ------------------------------
	Image impo = Image(im.width(),im.height(),im.channels());
	for (int x=0;x<im.width();x++){
		for (int y=0; y<im.height(); y++){
			for (int z=0; z<3; z++){
				impo(x,y,z)=1;
			}
		}
	}
	singleScalePaintImportance(im,impo,out,texture,size,N,noise);
	singleScalePaintImportance(im,sharpnessMap(im), out, texture,size/4,N,noise);
	return;
}
// Image gradientX(const Image &im, bool clamp) {
//   Filter sobelX(3, 3);
//   sobelX(0, 0) = -1.0;
//   sobelX(1, 0) = 0.0;
//   sobelX(2, 0) = 1.0;
//   sobelX(0, 1) = -2.0;
//   sobelX(1, 1) = 0.0;
//   sobelX(2, 1) = 2.0;
//   sobelX(0, 2) = -1.0;
//   sobelX(1, 2) = 0.0;
//   sobelX(2, 2) = 1.0;

//   Image imSobelX = sobelX.convolve(im, clamp);
//   return imSobelX;
// }

// Image gradientY(const Image &im, bool clamp) {

//   // sobel filtering in y direction
//   Filter sobelY(3, 3);
//   sobelY(0, 0) = -1.0;
//   sobelY(1, 0) = -2.0;
//   sobelY(2, 0) = -1.0;
//   sobelY(0, 1) = 0.0;
//   sobelY(1, 1) = 0.0;
//   sobelY(2, 1) = 0.0;
//   sobelY(0, 2) = 1.0;
//   sobelY(1, 2) = 2.0;
//   sobelY(2, 2) = 1.0;

//   Image imSobelY = sobelY.convolve(im, clamp);
//   return imSobelY;
// }

// Image computeTensor(const Image &im, float sigmaG, float factorSigma) {
//   // // --------- HANDOUT  PS07 ------------------------------
//   // Compute xx/xy/yy Tensor of an image. (stored in that order)
//   vector<Image> lumChrom = lumiChromi(im);
//   Image imb = gaussianBlur_separable(im,sigmaG);
//   Image imx = gradientX(imb,true);
//   Image imy = gradientY(imb,true);
//   Image res = Image(im.extent(0),im.extent(1),im.extent(2));
//   for (int x =0; x<res.width();x++){
//     for (int y=0; y<res.height();y++){
//       res(x,y,0) = powf(imx(x,y),2);
//       res(x,y,1) = imx(x,y)*imy(x,y);
//       res(x,y,2) = imy(x,y)*imy(x,y);
//     }
//   }
//   Image result = gaussianBlur_separable(res,sigmaG*factorSigma);

//   return result;
// }

// Image computeTensor(const Image &im, float sigmaG, float factorSigma) {
//  	// Compute xx/xy/yy Tensor of an image. (stored in that order)
//  	// // --------- HANDOUT  PS07 ------------------------------
//  	return Image(1,1,1);
// }


Image testAngle(const Image &im, float sigmaG, float factor) {
	// Extracts orientation of features in im. Angles should be mapped
	// to [0,1]
	// // --------- HANDOUT  PS12 ------------------------------
	Image tens = computeTensor(im,sigmaG,factor);
	Matrix M = Matrix::Zero(2,2);
	float ang;
	Eigen::EigenSolver<Matrix> eig;
	Image res = Image(im.width(),im.height(), im.channels());
	for (int x =0; x<res.width();x++){
		for (int y=0; y<res.height();y++){
			M(0,0) = tens(x,y,0);
			M(1,1) = tens(x,y,2);
			M(0,1) = tens(x,y,1);
			M(1,0) = tens(x,y,1);
			eig.compute(M);
			if (eig.eigenvalues()(0).real()<eig.eigenvalues()(1).real()){
				ang = atan2(-eig.eigenvectors().col(0)(1).real(),eig.eigenvectors().col(0)(0).real()); 
			}else{
				ang = atan2(-eig.eigenvectors().col(1)(1).real(),eig.eigenvectors().col(1)(0).real());
			}
			ang = ang<=0 ? ang+2*M_PI : ang;
			ang = ang/(2*M_PI);
			for (int z=0;z<3;z++){
				res(x,y,z) = ang;
			}
		}
	}

	return res;

}

vector<Image> rotateBrushes(const Image &im, int nAngles) {
	// helper function
	// Returns list of nAngles images rotated by 1*2pi/nAngles
	// // --------- HANDOUT  PS12 ------------------------------
	vector<Image> res;
	for (int i=0;i<nAngles;i++){
		res.push_back(rotate(im,2*M_PI*i/nAngles));
	}
	return res;
}

void singleScaleOrientedPaint(const Image &im, const Image &importance,
		Image &out, const Image &tensor, const Image &texture,int size, int N, 
		float noise, int nAngles) {
	// Similar to singleScalePaintImportant but brush strokes are oriented
	// according to tensor
	// // --------- HANDOUT  PS12 ------------------------------
	float factor = min(float(size)/texture.width(),float(size)/texture.height());
	Image tscaled = texture;
	float nfact = 0;
	int x, y, angN;
	Matrix M = Matrix::Zero(2,2);
	float ang;
	Eigen::EigenSolver<Matrix> eig;
	vector<float> color = {0,0,0};
	if (factor<1){
		tscaled = scaleLin(texture, factor);
	}
	vector<Image> rbrush = rotateBrushes(tscaled,nAngles);
	float sum=0;
	for (int x3=0;x3<importance.width();x3++){
		for (int y3=0;y3<importance.height();y3++){
			sum = sum + importance(x3,y3);
		}
	}
	int N1 = N/(sum/(importance.width()*importance.height()));


	for (int i=0;i<N1;i++){
		x = randm()*(out.width()-1);
		y = randm()*(out.height()-1);
		for (int z=0;z<3;z++){
			nfact = 1 - noise/2 + noise*randm();
			color[z] = nfact*im(x,y,z);
		}
		M(0,0) = tensor(x,y,0);
		M(1,1) = tensor(x,y,2);
		M(0,1) = tensor(x,y,1);
		M(1,0) = tensor(x,y,1);
		eig.compute(M);
		if (eig.eigenvalues()(0).real()<eig.eigenvalues()(1).real()){
			ang = atan2(-eig.eigenvectors().col(0)(1).real(),eig.eigenvectors().col(0)(0).real()); 
		}else{
			ang = atan2(-eig.eigenvectors().col(1)(1).real(),eig.eigenvectors().col(1)(0).real());
		}
		ang = ang<=0 ? ang+2*M_PI : ang;
		ang = ang/(2*M_PI);
		angN = floor(ang*nAngles);
		if (importance(x,y)>=randm()){
			brush(out,x,y,color,rbrush[angN]);
		}
	}

	return;
}

void orientedPaint(const Image &im, Image &out, const Image &texture, int N, int size, float noise) {
	// Similar to painterly() but strokes are oriented along the directions of maximal structure
	// // --------- HANDOUT  PS12 ------------------------------
	Image tensor = computeTensor(im);
	Image impo = Image(im.width(),im.height(),im.channels());
	for (int x=0;x<impo.width();x++){
		for (int y=0; y<impo.height(); y++){
			for (int z=0; z<3; z++){
				impo(x,y,z)=1;
			}
		}
	}
	singleScaleOrientedPaint(im,impo,out,tensor,texture,size,N,noise);
	singleScaleOrientedPaint(im,sharpnessMap(im),out,tensor,texture,size/4,N,noise);
	return;

}



