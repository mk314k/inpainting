/* -----------------------------------------------------------------
 * File:    Image.cpp
 * Created: 2015-08-29
 * Updated: 2019-08-10
 * -----------------------------------------------------------------
 *
 * The 6.815/6.865 Image Class
 *
 * ---------------------------------------------------------------*/

#include "Image.h"

using namespace std;

// --------- HANDOUT  PS04 ------------------------------
// ------------------------------------------------------

// obtain minimum pixel value
float Image::min() const {
  float minf = FLT_MAX;
  for (int i = 0; i < number_of_elements(); i++) {
    minf = std::min(minf, (*this)(i));
  }
  return minf;
}

// obtain maximum pixel value
float Image::max() const {
  float maxf = -FLT_MAX;
  for (int i = 0; i < number_of_elements(); i++) {
    maxf = std::max(maxf, (*this)(i));
  }
  return maxf;
}
// --------- END of PS04 --------------------------------

// --------- HANDOUT  PS02 ------------------------------
// ------------------------------------------------------

// Safe Accessor that will return a black pixel (clamp = false) or the nearest
// pixel value (clamp = true) when indexing out of the bounds of the image
float Image::smartAccessor(int x, int y, int z, bool clamp) const {
  // // --------- HANDOUT  PS02 ------------------------------
  // return 0.0f; // change this

  // --------- SOLUTION PS02 ------------------------------
  float black = 0.0;
  int x0 = x;
  int y0 = y;

  if (y >= height()) {
    if (clamp) {
      y0 = height() - 1;
    } else {
      return black;
    }
  }

  if (y < 0) {
    if (clamp) {
      y0 = 0;
    } else {
      return black;
    }
  }

  if (x >= width()) {
    if (clamp) {
      x0 = width() - 1;
    } else {
      return black;
    }
  }

  if (x < 0) {
    if (clamp) {
      x0 = 0;
    } else {
      return black;
    }
  }

  return (*this)(x0, y0, z);
}
// ---------------- END of PS02 -------------------------------------

// --------- HANDOUT  PS01 ------------------------------

long long Image::number_of_elements() const {
  // --------- HANDOUT  PS01 ------------------------------
  // returns the number of elements in the image.
  // An RGB (3 color channels) image of 100 × 100 pixels has 30000 elements
  // return 0; // change this

  // --------- SOLUTION PS01 ------------------------------
  return image_data.size();
}

// -------------- Accessors and Setters -----------------
const float &Image::operator()(int x) const {
  // --------- HANDOUT  PS01 ------------------------------
  // Linear accessor to the image data
  // throw NotImplementedException(); // change this

  // --------- SOLUTION PS01 ------------------------------
  if (x < 0 || x >= number_of_elements())
    throw OutOfBoundsException();
  return image_data[x];
}

const float &Image::operator()(int x, int y) const {
  // --------- HANDOUT  PS01 ------------------------------
  // Accessor to the image data at channel 0
  // throw NotImplementedException(); // change this

  // --------- SOLUTION PS01 ------------------------------
  if (x < 0 || x >= width())
    throw OutOfBoundsException();
  if (y < 0 || y >= height())
    throw OutOfBoundsException();

  return image_data[x * stride_[0] + y * stride_[1]];
}

const float &Image::operator()(int x, int y, int z) const {
  // --------- HANDOUT  PS01 ------------------------------
  // Accessor to the image data at channel z
  // throw NotImplementedException(); // change this

  // --------- SOLUTION PS01 ------------------------------
  if (x < 0 || x >= width())
    throw OutOfBoundsException();
  if (y < 0 || y >= height())
    throw OutOfBoundsException();
  if (z < 0 || z >= channels())
    throw OutOfBoundsException();

  return image_data[x * stride_[0] + y * stride_[1] + stride_[2] * z];
}

float &Image::operator()(int x) {
  // --------- HANDOUT  PS01 ------------------------------
  // Linear setter to the image data
  // throw NotImplementedException(); // change this

  // --------- SOLUTION PS01 ------------------------------
  if (x < 0 || x >= number_of_elements())
    throw OutOfBoundsException();
  return image_data[x];
}

float &Image::operator()(int x, int y) {
  // --------- HANDOUT  PS01 ------------------------------
  // Setter to the image data at channel 0
  // throw NotImplementedException(); // change this

  // --------- SOLUTION PS01 ------------------------------
  if (x < 0 || x >= width())
    throw OutOfBoundsException();
  if (y < 0 || y >= height())
    throw OutOfBoundsException();

  return image_data[x * stride_[0] + y * stride_[1]];
}

float &Image::operator()(int x, int y, int z) {
  // --------- HANDOUT  PS01 ------------------------------
  // Setter to the image data at channel z
  // throw NotImplementedException(); // change this

  // --------- SOLUTION PS01 ------------------------------
  if (x < 0 || x >= width())
    throw OutOfBoundsException();
  if (y < 0 || y >= height())
    throw OutOfBoundsException();
  if (z < 0 || z >= channels())
    throw OutOfBoundsException();

  return image_data[x * stride_[0] + y * stride_[1] + stride_[2] * z];
}

void Image::set_color(float r, float g, float b) {
  // --------- HANDOUT  PS01 ------------------------------
  // Set the image pixels to the corresponding values
  // throw NotImplementedException(); // change this

  // --------- SOLUTION PS01 ------------------------------
  for (int i = 0; i < width() * height(); ++i) {
    image_data[i] = r;
    if (channels() > 1) // have second channel
      image_data[i + stride_[2]] = g;
    if (channels() > 2) // have third channel
      image_data[i + 2 * stride_[2]] = b;
  }
}

void Image::create_rectangle(int xstart, int ystart, int xend, int yend,
                             float r, float g, float b) {
  // --------- HANDOUT  PS01 ------------------------------
  // Set the pixels inside the rectangle to the specified color
  // throw NotImplementedException(); // change this

  // --------- SOLUTION PS01 ------------------------------
  if (xstart < 0 || xstart >= width() || ystart < 0 || ystart >= height())
    throw OutOfBoundsException();
  if (xend < 0 || xend >= width() || yend < 0 || yend >= height())
    throw OutOfBoundsException();

  float col[3] = {r, g, b};
  for (int w = xstart; w <= xend; ++w) {
    for (int h = ystart; h <= yend; ++h) {
      for (int c = 0; c < channels(); ++c) {
        (*this)(w, h, c) = col[c];
      }
    }
  }
}

void Image::create_line(int xstart, int ystart, int xend, int yend, float r,
                        float g, float b) {
  // --------- HANDOUT  PS01 ------------------------------
  // Create a line segment with specified color
  // throw NotImplementedException(); // change this

  // --------- SOLUTION PS01 ------------------------------
  if (xstart < 0 || xstart >= width() || ystart < 0 || ystart >= height())
    throw OutOfBoundsException();
  if (xend < 0 || xend >= width() || yend < 0 || yend >= height())
    throw OutOfBoundsException();

  int valid_channels = channels() > 3 ? 3 : channels();
  float col[3] = {r, g, b};
  float x = xstart;
  float y = ystart;
  int delta_x = xend - xstart;
  int delta_y = yend - ystart;
  int delta = std::max(std::abs(delta_x), std::abs(delta_y));
  for (int i = 0; i <= delta; i++) {
    int ix = std::round(x), iy = std::round(y);
    if (ix >= 0 && ix < width() && iy >= 0 && iy < height()) {
      for (int c = 0; c < valid_channels; ++c) {
        (*this)(ix, iy, c) = col[c];
      }
    }
    x += (float(delta_x) / float(delta));
    y += (float(delta_y) / float(delta));
  }
}

// ---------------- END of PS01 -------------------------------------

/*********************************************************************
 *                    DO NOT EDIT BELOW THIS LINE                    *
 *********************************************************************/

int Image::debugWriteNumber = 0;

Image::Image(int width, int height, int channels, const std::string &name_) {
  initialize_image_metadata(width, height, channels, name_);
  // Initialize image data
  long long size_of_data = 1;
  for (unsigned int k = 0; k < DIMS; k++) {
    size_of_data *= dim_values[k];
  }
  image_data = std::vector<float>(size_of_data, 0.f);
}

void Image::initialize_image_metadata(int w, int h, int c,
                                      const std::string &name_) {
  if (w < 1)
    throw IllegalDimensionException();
  if (h < 1)
    throw IllegalDimensionException();
  if (c != 1 && c != 3)
    throw IllegalDimensionException();

  dim_values[0] = w;
  dim_values[1] = h;
  dim_values[2] = c;
  stride_[0] = 1;
  stride_[1] = w;
  stride_[2] = w * h;
  image_name = name_;
}

Image::Image(const std::string &filename) {
  std::vector<unsigned char> uint8_image;
  unsigned int h;
  unsigned int w;
  unsigned int c = 4;
  unsigned int oc = 3; // Throw away transparency

  // Decode PNG file In column major order with packed color values
  unsigned err = lodepng::decode(uint8_image, w, h, filename.c_str());
  if (err == 48) {
    throw FileNotFoundException();
  }

  image_data = std::vector<float>(h * w * oc, 0);

  for (unsigned int x = 0; x < w; x++) {
    for (unsigned int y = 0; y < h; y++) {
      for (unsigned int z = 0; z < oc; z++) {
        image_data[x + y * w + z * w * h] =
            uint8_to_float(uint8_image[z + x * c + y * c * w]);
      }
    }
  }

  initialize_image_metadata(w, h, oc, filename);
}

Image::~Image() {} // Nothing to clean up

void Image::write(const std::string &filename) const {
  if (channels() != 1 && channels() != 3 && channels() != 4)
    throw ChannelException();
  int png_channels = 4;
  std::vector<unsigned char> uint8_image(height() * width() * png_channels,
                                         255);
  int z;
  for (int x = 0; x < width(); x++) {
    for (int y = 0; y < height(); y++) {
      for (z = 0; z < channels(); z++) {
        uint8_image[z + x * png_channels + y * png_channels * width()] =
            float_to_uint8(
                image_data[x + y * width() + z * width() * height()]);
      }
      for (; z < 3; z++) { // Only executes when there is one channel
        uint8_image[z + x * png_channels + y * png_channels * width()] =
            float_to_uint8(
                image_data[x + y * width() + 0 * width() * height()]);
      }
    }
  }
  lodepng::encode(filename.c_str(), uint8_image, width(), height());
}

void Image::debug_write() const {
  std::ostringstream ss;
  ss << "./Output/" << debugWriteNumber << ".png";
  std::string filename = ss.str();
  write(filename);
  debugWriteNumber++;
}

float Image::uint8_to_float(const unsigned char &in) {
  return ((float)in) / (255.0f);
}

unsigned char Image::float_to_uint8(const float &in) {
  float out = in;
  if (out < 0)
    out = 0;
  if (out > 1)
    out = 1;
  return (unsigned char)(255.0f * out);
}

void compareDimensions(const Image &im1, const Image &im2) {
  for (int i = 0; i < 3; i++) {
    if (im1.extent(i) != im2.extent(i))
      throw MismatchedDimensionsException();
  }
}

Image operator+(const Image &im1, const Image &im2) {
  compareDimensions(im1, im2);
  long long total_pixels = im1.number_of_elements();

  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    output(i) = im1(i) + im2(i);
  }
  return output;
}

Image operator-(const Image &im1, const Image &im2) {
  compareDimensions(im1, im2);
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    output(i) = im1(i) - im2(i);
  }
  return output;
}

Image operator*(const Image &im1, const Image &im2) {
  compareDimensions(im1, im2);
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    output(i) = im1(i) * im2(i);
  }
  return output;
}

Image operator/(const Image &im1, const Image &im2) {
  compareDimensions(im1, im2);
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    if (im2(i) == 0)
      throw DivideByZeroException();
    output(i) = im1(i) / im2(i);
  }
  return output;
}

Image operator+(const Image &im1, const float &c) {
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    output(i) = im1(i) + c;
  }
  return output;
}

Image operator-(const Image &im1, const float &c) {
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    output(i) = im1(i) - c;
  }
  return output;
}
Image operator*(const Image &im1, const float &c) {
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    output(i) = im1(i) * c;
  }
  return output;
}
Image operator/(const Image &im1, const float &c) {
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  if (c == 0)
    throw DivideByZeroException();
  for (int i = 0; i < total_pixels; i++) {
    output(i) = im1(i) / c;
  }
  return output;
}

Image operator+(const float &c, const Image &im1) {
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    output(i) = im1(i) + c;
  }
  return output;
}

Image operator-(const float &c, const Image &im1) {
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    output(i) = c - im1(i);
  }
  return output;
}

Image operator*(const float &c, const Image &im1) {
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    output(i) = im1(i) * c;
  }
  return output;
}
Image operator/(const float &c, const Image &im1) {
  long long total_pixels = im1.number_of_elements();
  Image output(im1.extent(0), im1.extent(1), im1.extent(2));
  for (int i = 0; i < total_pixels; i++) {
    if (im1(i) == 0)
      throw DivideByZeroException();
    output(i) = c / im1(i);
  }
  return output;
}
