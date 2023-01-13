/* --------------------------------------------------------------------------
 * File:    a12_main.cpp
 * Created: 2021-10-28
 * --------------------------------------------------------------------------
 *
 *
 *
 * ------------------------------------------------------------------------*/

#include "Image.h"
#include "a9.h"
#include "basicImageManipulation.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "matrix.h"

using namespace std;

Image maskThreshold(const Image &mask){
  Image res = Image(mask.width(),mask.height(),3);
  float val =0.0;
  for (int x=0;x<mask.width();x++){
      for (int y=0; y<mask.height();y++){
        if (mask(x,y)>=0.5){
          val=1.0;
        }else{
          val = 0.0;
        }
        for (int z=0;z<3;z++){
            res(x,y,z) = val;
        }
      }
  }
  return res;
}
void imageCopy(const Image &s, Image &d, int y1, int x1, int y2=0, int x2=0){
  int w = min(d.width(),s.width());
  int h = min(d.height(),s.height());
  for (int x=0;x<w;x++){
    for (int y=0; y<h;y++){
      for (int z=0;z<d.channels();z++){
        d(x+x2,y+y2,z) = s(x+x1,y+y1,z);
      }
    }
  }
}
void testDeconvCG_reg(){
  Image fg("./Input/bear.png");
  vector<float> fData = gauss2DFilterValues(3.0, 1.0);
  int k = sqrt(fData.size());
  Filter gauss(fData, k, k);
  Image blur = gauss.convolve(fg);
  blur.write("./Output/bear_blur.png");
  Image out = deconvCG_reg(blur,gauss);
  out.write("./Output/deconv_cg_reg.png");
}
void testPoisson(){
  Image fg("./Input/bear.png");
  Image bg("./Input/waterpool.png");
  Image mask("./Input/mask.png");
  mask = maskThreshold(mask);
  Image bgn = Image(fg.width(),fg.height(),fg.channels());
  imageCopy(bg,bgn,50,10);
  Image out = poisson(bgn,fg,mask,200);
  out.write("./Output/poisson.png");
}

void testPoissonCG(){
  Image fg("./Input/bear.png");
  Image bg("./Input/waterpool.png");
  Image mask("./Input/mask.png");
  mask = maskThreshold(mask);
  Image bgn = Image(fg.width(),fg.height(),fg.channels());
  imageCopy(bg,bgn,50,10);
  Image out = poissonCG(bgn,fg,mask,100);
  imageCopy(out,bg,0,0,50,10);
  bg.write("./Output/poisson_cg.png");
}
void testFlatPoissonCG(){
  Image bg("./Input/waterpool.png");
  Image mask("./Input/mask.png");
  Image fg = Image(mask.width(), mask.height(), bg.channels());
  mask = maskThreshold(mask);
  Image bgn = Image(fg.width(),fg.height(),fg.channels());
  imageCopy(bg,bgn,50,10);
  Image out = poissonCG(bgn,fg,mask,200);
  out.write("./Output/poisson_flat_cg.png");
}
void testLaplacian(){
  Image fg("./Input/bear.png");
  Image out = applyLapacian(fg);
  out.write("./Output/bear_laplace.png");
}
void testInpaint(){
  Image bg("./Input/bhanju.png");
  Image mask = Image(bg.width(), bg.height(),bg.channels());
  bg.create_line(10,10,50,50,0.4,0.7,0.2);
  mask.create_line(10,10,50,50,1.0,1.0,1.0);
  bg.create_rectangle(60,60,70,200,0.4,0.7,0.2);
  mask.create_rectangle(60,60,70,200,1.0,1.0,1.0);
  mask.write("./Output/mask.png");
  bg.write("./Output/masked.png");
  Image out = inpaint_poisson(bg,mask);
  out.write("./Output/inpaint_cg.png");
}
void testTensor(){
  Image bg("./Input/bhanju.png");
  Image out = computeTensor(bg);
  out.write("./Output/tensor.png");
}
// This is a way for you to test your functions.
// We will only grade the contents of npr.cpp
int main() {
  cout<<"nothing done in a9_main.cpp, debug me !"<< endl;
  // testPoisson();
  testPoissonCG();
  // testFlatPoissonCG();
  // testDeconvCG_reg();
  // testLaplacian();
  // testInpaint();
  // testTensor();
  return 0;
}
