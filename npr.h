/* --------------------------------------------------------------------------
 * File:    npr.h
 * Created: 2021-10-28
 * --------------------------------------------------------------------------
 *
 *  Assignment 8
 *
 * ------------------------------------------------------------------------*/

#ifndef __NPR__H
#define __NPR__H

#include "Image.h"
#include "basicImageManipulation.h"
#include <iostream>
#include <math.h>

using namespace std;

void brush(Image &im, int x, int y, vector<float> color, const Image &texture);
void singleScalePaint(const Image &im, Image &out, const Image &texture, int size=10, int N=10000, float noise=0.3);
void singleScalePaintImportance(const Image &im, const Image &importance, Image &out, const Image &texture, 
			int size=10, int N=10000, float noise=0.3);

Image sharpnessMap(const Image &im, float sigma=1.0f);
void painterly(const Image &im, Image &out, const Image &texture, int N=10000, int size=50, float noise=0.3);

Image computeTensor(const Image &im, float sigmaG=1.0f, float factor=4.0f);
Image getBlurredLumi(const Image &im, float sigmaG);

Image testAngle(const Image &im, float sigmaG=3.0f, float factor=5.0f);

vector<Image> rotateBrushes(const Image &im, int nAngles);

void singleScaleOrientedPaint(const Image &im, const Image &importance, Image &out, const Image &tensor,
	const Image &texture, int size, int N, float noise, int nAngles=36);

void orientedPaint(const Image &im, Image &out, const Image &texture, int N=10000, int size=50, float noise=0.3);


#endif
