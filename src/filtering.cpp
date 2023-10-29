

#include "filtering.h"
#include <cassert>
#include <cmath>

using namespace std;

Image boxBlur(const Image &im, int k, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Convolve an image with a box filter of size k by k
  // It is safe to asssume k is odd.
  // return im; // change this

  // --------- SOLUTION PS02 ------------------------------
  // Create a new empty image
  Image filtered(im.width(), im.height(), im.channels());
  int sideSize = int((k - 1.0f) / 2.0f);
  float normalizer = 1.0f / float(k * k);
  float accum = 0.0f;

  // For every pixel in the image
  for (int z = 0; z < filtered.channels(); z++) {
    for (int y = 0; y < filtered.height(); y++) {
      for (int x = 0; x < filtered.width(); x++) {
        // Accumulate the sum in the pixel's kxk neighborhood
        accum = 0.0f;
        for (int yBox = -sideSize; yBox < -sideSize + k; yBox++) {
          for (int xBox = -sideSize; xBox < -sideSize + k; xBox++) {
            accum += im.smartAccessor(x - xBox, y - yBox, z, clamp);
          }
        }

        // Assign normalized result of convolution to output pixel
        filtered(x, y, z) = accum * normalizer;
      }
    }
  }
  return filtered;
}

Image Filter::convolve(const Image &im, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Write a convolution function for the filter class
  // return im; // change this

  // --------- SOLUTION PS02 ------------------------------
  Image imFilter(im.width(), im.height(), im.channels());

  int sideW = int((width - 1.0) / 2.0);
  int sideH = int((height - 1.0) / 2.0);
  float accum;

  // For every pixel in the image
  for (int z = 0; z < imFilter.channels(); z++) {
    for (int y = 0; y < imFilter.height(); y++) {
      for (int x = 0; x < imFilter.width(); x++) {
        accum = 0.0;
        for (int yFilter = 0; yFilter < height; yFilter++) {
          for (int xFilter = 0; xFilter < width; xFilter++) {
            // Sum the image pixel values weighted by the filter
            // flipped kernel, xFilter, yFilter have different
            // signs in filter and im
            accum += operator()(xFilter, yFilter) *
                     im.smartAccessor(x - xFilter + sideW, y - yFilter + sideH,
                                      z, clamp);
          }
        }
        // Assign the pixel the value from convolution
        imFilter(x, y, z) = accum;
      }
    }
  }
  return imFilter;
}
Filter Filter::transpose(){
    Filter res = Filter(height,width);
    for (int y=0; y<height; y++) {
        for (int x=0; x<width; x++) {
            res(y,x) = operator()(x,y);
        }
    }
    return res;
}

Image boxBlur_filterClass(const Image &im, int k, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Reimplement the box filter using the filter class.
  // check that your results match those in the previous function "boxBlur"
  // return im; // change this

  // --------- SOLUTION PS02 ------------------------------
  vector<float> fData(k * k, 1.0 / (k * k));
  Filter boxFilter(fData, k, k);
  Image imFilter = boxFilter.convolve(im, clamp);
  return imFilter;
}

Image gradientMagnitude(const Image &im, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Uses a Sobel kernel to compute the horizontal and vertical components
  // of the gradient of an image and returns the gradient magnitude.
  // return im; // change this

  // --------- SOLUTION PS02 ------------------------------
  // Sobel filtering in x direction
  Filter sobelX(vector<float>{-1.0, 0.0, 1.0, -2.0, 0.0, 2.0, -1.0, 0.0, 1.0},
                3, 3);

  Image imSobelX = sobelX.convolve(im, clamp);

  // Sobel filtering in y direction
  Filter sobelY(vector<float>{-1.0, -2.0, -1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0},
                3, 3);
  Image imSobelY = sobelY.convolve(im, clamp);

  // squared magnitude
  Image magnitude = imSobelX * imSobelX + imSobelY * imSobelY;

  // take the square root
  for (int i = 0; i < magnitude.number_of_elements(); i++) {
    magnitude(i) = sqrt(magnitude(i));
  }
  return magnitude;
}

vector<float> gauss1DFilterValues(float sigma, float truncate) {
  // --------- HANDOUT  PS02 ------------------------------
  // Create a vector containing the normalized values in a 1D Gaussian filter
  // Truncate the gaussian at truncate*sigma.
  // return vector<float>();

  // --------- SOLUTION PS02 ------------------------------
  // compute the support of the filter
  int offset = int(ceil(truncate * sigma));
  int filterSize = 2 * offset + 1;

  vector<float> fData(filterSize, 0.0f);

  // compute the un-normalized value of the gaussian
  float normalizer = 0.0f;
  for (int i = 0; i < filterSize; i++) {
    fData[i] = exp(-(i - offset) * (i - offset) / (2.0f * sigma * sigma));
    normalizer += fData[i];
  }

  // normalize
  for (int i = 0; i < filterSize; i++) {
    fData[i] /= normalizer;
  }
  return fData;
}

Image gaussianBlur_horizontal(const Image &im, float sigma, float truncate,
                              bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Gaussian blur across the rows of an image
  // return im;

  // --------- SOLUTION PS02 ------------------------------
  // Filter in the x direction
  vector<float> fData = gauss1DFilterValues(sigma, truncate);
  Filter gaussX(fData, fData.size(), 1);
  Image imFilter = gaussX.convolve(im, clamp);
  return imFilter;
}

vector<float> gauss2DFilterValues(float sigma, float truncate) {
  // --------- HANDOUT  PS02 ------------------------------
  // Create a vector containing the normalized values in a 2D Gaussian
  // filter. Truncate the gaussian at truncate*sigma.
  // return vector<float>();

  // --------- SOLUTION PS02 ------------------------------
  // Compute the filter size
  int offset = int(ceil(truncate * sigma));
  int k = 2 * offset + 1;
  int filterSize = k * k;

  vector<float> fData(filterSize, 0);

  int idx = 0;
  float normalizer = 0.0;

  // compute the unnormalized value of the gaussian and put it in a
  // row-major vector
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      fData[idx] =
          exp(-((i - offset) * (i - offset) + (j - offset) * (j - offset)) /
              (2.0 * sigma * sigma));
      normalizer += fData[idx];
      idx++;
    }
  }

  // normalize
  for (int i = 0; i < filterSize; i++) {
    fData[i] /= normalizer;
  }

  return fData;
}

Image gaussianBlur_2D(const Image &im, float sigma, float truncate,
                      bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Blur an image with a full 2D rotationally symmetric Gaussian kernel
  // return im;

  // --------- SOLUTION PS02 ------------------------------
  // Blur using a 2D gaussian filter
  vector<float> fData = gauss2DFilterValues(sigma, truncate);
  int k = sqrt(fData.size());
  Filter gauss(fData, k, k);
  Image imFilter = gauss.convolve(im, clamp);

  return imFilter;
}

Image gaussianBlur_separable(const Image &im, float sigma, float truncate,
                             bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Use principles of separability to blur an image using 2 1D Gaussian
  // Filters
  // return im;

  // --------- SOLUTION PS02 ------------------------------
  // Blur using two 1D filters in the x and y directions
  vector<float> fData = gauss1DFilterValues(sigma, truncate);
  Filter gaussX(fData, fData.size(), 1);
  Filter gaussY(fData, 1, fData.size());
  Image imFilter = gaussX.convolve(im, clamp);
  imFilter = gaussY.convolve(imFilter, clamp);

  return imFilter;
}

Image unsharpMask(const Image &im, float sigma, float truncate, float strength,
                  bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Sharpen an image
  // return im;

  // --------- SOLUTION PS02 ------------------------------
  // Get the low pass image
  Image lowPass = gaussianBlur_separable(im, sigma, truncate, clamp);
  // Subtract it from the original image to get the high pass image
  Image highPass = im - lowPass;

  // Increase the highPass component to sharpen
  Image sharp = im + strength * highPass;
  return sharp;
}

Image bilateral(const Image &im, float sigmaRange, float sigmaDomain,
                float truncateDomain, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Denoise an image using the bilateral filter
  // return im;

  // --------- SOLUTION PS02 ------------------------------
  Image imFilter(im.width(), im.height(), im.channels());

  // calculate the filter size
  int offset = int(ceil(truncateDomain * sigmaDomain));
  int sizeFilt = 2 * offset + 1;
  float accum, tmp, range_dist, normalizer, factorDomain, factorRange;

  // for every pixel in the image
  for (int z = 0; z < imFilter.channels(); z++)
    for (int y = 0; y < imFilter.height(); y++)
      for (int x = 0; x < imFilter.width(); x++) {
        // initilize normalizer and sum value to 0 for every pixel location
        normalizer = 0.0f;
        accum = 0.0f;

        // sum over the filter's support
        for (int yFilter = 0; yFilter < sizeFilt; yFilter++)
          for (int xFilter = 0; xFilter < sizeFilt; xFilter++) {
            // calculate the distance between the 2 pixels (in range)
            range_dist = 0.0f; // |R-R1|^2 + |G-G1|^2 + |B-B1|^2
            for (int z1 = 0; z1 < imFilter.channels(); z1++) {
              tmp = im.smartAccessor(x, y, z1, clamp); // center pixel
              tmp -=
                  im.smartAccessor(x + xFilter - offset, y + yFilter - offset,
                                   z1, clamp); // neighbor
              tmp *= tmp;                      // square
              range_dist += tmp;
            }

            // calculate the exponential weight from the domain and range
            factorDomain = exp(-((xFilter - offset) * (xFilter - offset) +
                                 (yFilter - offset) * (yFilter - offset)) /
                               (2.0 * sigmaDomain * sigmaDomain));
            factorRange = exp(-range_dist / (2.0 * sigmaRange * sigmaRange));

            normalizer += factorDomain * factorRange;
            accum += factorDomain * factorRange *
                     im.smartAccessor(x + xFilter - offset,
                                      y + yFilter - offset, z, clamp);
          }

        // set pixel in filtered image to weighted sum of values in the
        // filter region
        imFilter(x, y, z) = accum / normalizer;
      }

  return imFilter;
}

Image bilaYUV(const Image &im, float sigmaRange, float sigmaY, float sigmaUV,
              float truncateDomain, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // 6.865 only
  // Bilateral Filter an image separately for
  // the Y and UV components of an image
  // return im;

  // --------- SOLUTION PS02 ------------------------------
  // convert from RGB to YUV
  Image imYUV = rgb2yuv(im);

  // We pass the whole imYUV to bilateral, since we want to compute the
  // weight on the full YUV range
  Image bilY = bilateral(imYUV, sigmaRange, sigmaY, truncateDomain, clamp);
  Image bilUV = bilateral(imYUV, sigmaRange, sigmaUV, truncateDomain, clamp);

  // put the Y and UV parts of the image back into one image
  for (int i = 0; i < im.width(); i++) {
    for (int j = 0; j < im.height(); j++) {
      imYUV(i, j, 0) = bilY(i, j, 0);
      imYUV(i, j, 1) = bilUV(i, j, 1);
      imYUV(i, j, 2) = bilUV(i, j, 2);
    }
  }

  // convert from YUV back to RGB
  Image bilRGB = yuv2rgb(imYUV);
  return bilRGB;
}

/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// Create an image of 0's with a value of 1 in the middle. This function
// can be used to test that you have properly set the kernel values in your
// Filter object. Make sure to set k to be larger than the size of your kernel
Image impulseImg(int k) {
  // initlize a kxkx1 image of all 0's
  Image impulse(k, k, 1);

  // set the center pixel to have intensity 1
  int center = floor(k / 2);
  impulse(center, center, 0) = 1.0f;

  return impulse;
}

// ------------- FILTER CLASS -----------------------
Filter::Filter(const vector<float> &fData, int fWidth, int fHeight)
    : kernel(fData), width(fWidth), height(fHeight) {
  assert(fWidth * fHeight == (int)fData.size());
}

Filter::Filter(int fWidth, int fHeight)
    : kernel(std::vector<float>(fWidth * fHeight, 0)), width(fWidth),
      height(fHeight) {}

Filter::~Filter() {}

const float &Filter::operator()(int x, int y) const {
  if (x < 0 || x >= width)
    throw OutOfBoundsException();
  if (y < 0 || y >= height)
    throw OutOfBoundsException();

  return kernel[x + y * width];
}

float &Filter::operator()(int x, int y) {
  if (x < 0 || x >= width)
    throw OutOfBoundsException();
  if (y < 0 || y >= height)
    throw OutOfBoundsException();

  return kernel[x + y * width];
}
// --------- END FILTER CLASS -----------------------
Image gradientX(const Image &im, bool clamp){
    Filter sobelX(3, 3);
    sobelX(0,0) = -1.0; sobelX(1,0) = 0.0; sobelX(2,0) = 1.0;
    sobelX(0,1) = -2.0; sobelX(1,1) = 0.0; sobelX(2,1) = 2.0;
    sobelX(0,2) = -1.0; sobelX(1,2) = 0.0; sobelX(2,2) = 1.0;

    Image imSobelX = sobelX.convolve(im, clamp);
    return imSobelX;
}


Image gradientY(const Image &im, bool clamp) {

    // sobel filtering in y direction
    Filter sobelY(3, 3);
    sobelY(0,0) = -1.0; sobelY(1,0) = -2.0; sobelY(2,0) = -1.0;
    sobelY(0,1) = 0.0; sobelY(1,1) = 0.0; sobelY(2,1) = 0.0;
    sobelY(0,2) = 1.0; sobelY(1,2) = 2.0; sobelY(2,2) = 1.0;


    Image imSobelY = sobelY.convolve(im, clamp);
    return imSobelY;
}

Image maximum_filter(const Image &im, float maxiDiam) {
    float mi = floor((maxiDiam) / 2);
    float ma = maxiDiam - mi - 1;

    Image mf(im.width(), im.height(), im.channels());
    for (int c = 0; c < im.channels(); c++)
    for (int i = mi; i < im.width() - ma; i++)
    for (int j = mi; j < im.height() - ma; j++)
    {

        for (int deli = -mi; deli <= ma; deli++)
        for (int delj = -mi; delj <= ma; delj++)
        {
            mf(i, j, c) = max(mf(i, j, c), im(i+deli, j+delj, c));
        }

    }
    return mf;
}
