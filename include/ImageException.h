/* -----------------------------------------------------------------
 * File:    ImageException.h
 * Created: 2015-08-29
 * Updated: 2019-08-10
 * -----------------------------------------------------------------
 *
 * 6.815/6.865 Image Exceptions
 *
 * ---------------------------------------------------------------*/

#ifndef __IMAGEEXCEPTION__H
#define __IMAGEEXCEPTION__H

#include <stdexcept>

class DivideByZeroException : public std::runtime_error {
public:
  DivideByZeroException() : std::runtime_error("Divisor is zero") {}
};

class MismatchedDimensionsException : public std::runtime_error {
public:
  MismatchedDimensionsException()
      : std::runtime_error("Image dimensions are not the same.") {}
};

class IllegalDimensionException : public std::runtime_error {
public:
  IllegalDimensionException()
      : std::runtime_error("Image dimensions must be valid.") {}
};

class ChannelException : public std::runtime_error {
public:
  ChannelException()
      : std::runtime_error("Number of channels must be 1, 3 or 4 when \
                    writing to image.") {}
};

class OutOfBoundsException : public std::runtime_error {
public:
  OutOfBoundsException()
      : std::runtime_error("Index is out of the image bounds.") {}
};

class InvalidArgument : public std::runtime_error {
public:
  InvalidArgument() : std::runtime_error("Argument is not valid.") {}
};

class FileNotFoundException : public std::runtime_error {
public:
  FileNotFoundException()
      : std::runtime_error("Empty input or file does not exist.") {}
};

class NotImplementedException : public std::runtime_error {
public:
  NotImplementedException() : std::runtime_error("Method not implemented") {}
};

#endif
