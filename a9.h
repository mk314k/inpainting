#ifndef __A9_H
#define __A9_H

#include "Image.h"
#include "basicImageManipulation.h"
#include <iostream>
#include <math.h>
#include "filtering.h"

using namespace std;


void someFunction();
float dotIm(const Image &im1, const Image &im2, const bool &zeroClamp = true);
Image applyKernel(const Image &im, Filter &filt);
Image applyConjugatedKernel(const Image &im, Filter &filt);
Image computeResidual(const Image &x, const Image &y, Filter &kernel);
Image applyRegularizedOperator(const Image &x, Filter &kernel, const float &lamb);
Image applyLapacian(const Image &im);


Image deconvCG_reg(const Image &y, Filter &kernel, const float &lamb= 0.05, const int &niter=10);

Image poisson(const Image &bg, const Image &fg, const Image &mask, const int niter);
Image poissonCG(const Image &bg, const Image &fg, const Image &mask, const int niter=200);
Image inpaint_poisson(const Image &im, const Image &mask);

vector<Image> channel_split(const Image &im);
Image channel_combine(const vector<Image> &IM);
Image channel_copy(const Image &im);

Image computeTensor(const Image &im, float sigmaG=1.0f, float factorSigma=4.0f) ;

float solve(int i1,int j1,int i2,int j2, Image &f, Image &T,int KNOWN = 0);

void inpaint_point(int i,int j, const Image &mask, const Image &T, const Image &Tx, const Image &Ty, Image &out, int eps=2);
Image inpaint(const Image &im, const Image &mask);

#endif /* end of include guard: A10_H_PHUDVTKB */

