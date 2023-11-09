
[![C/C++ CI](https://github.com/mk314k/ImgLib/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/mk314k/ImgLib/actions/workflows/c-cpp.yml)
![License](https://img.shields.io/github/license/mk314k/ImgLib)
![Version](https://img.shields.io/badge/version-1.0.0-blue)

# Inpainting 
ImgLib is a comprehensive image processing library developed from scratch, providing a wide range of features and capabilities for image manipulation. Whether you're working on image analysis, computer vision, or graphics applications, ImgLib is a valuable tool with the following features:

Png Read/Write: ImgLib supports reading and writing PNG images. It leverages the power of LodePNG (version 20141130, Copyright (c) 2005-2014 Lode Vandevenne) for these operations.

Arithmetic Operators on Images: Perform mathematical operations on images to create stunning visual effects.

Brightness and Contrast Adjustments: Easily adjust the brightness and contrast of images for better visualization.

Luminance and Chrominance Extraction: Separate the luminance and chrominance components of an image, which is useful for various image processing tasks.

Color Spaces Conversion: ImgLib provides conversion between different color spaces, including YUV and RGB.

Gamma Coding: Apply gamma correction to images to enhance or modify their visual appearance.

Image Scaling: Resize images using various algorithms such as Nearest Neighbor, Linear, Bicubic, and Lanczos for different scaling needs.

Image Rotation: Rotate images at different angles to achieve the desired orientation.

Image Filtering: Apply various image filters, including convolution, blurring (Box Blur, Gaussian Blur), sharpening, bilateral filtering, gradients, and MaxPool operations.

Image Texturing and Brush Painting: Add textures to images and create artistic effects with brush painting.

Image Morphing: Achieve smooth transitions between images, which is valuable for animations and visual effects.

## Dependencies
ImgLib includes the LodePNG library (Copyright (c) 2005-2014 Lode Vandevenne) as part of the project, enabling PNG image read and write functionalities.

## Ongoing Work
Our ongoing development efforts are focused on creating CUDA versions of the features to harness the power of GPU acceleration for faster image processing.

## Getting Started
To run ImgLib locally, follow these simple steps:

Clone this Repository: Begin by cloning this repository to your local machine.

Compile and Run: Use the provided build system to compile and run ImgLib. For example, you can use the command make run to start using the library.
