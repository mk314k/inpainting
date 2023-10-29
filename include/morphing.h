#ifndef __MORPHING__H
#define __MORPHING__H

#include "Image.h"
#include "basicImageManipulation.h"
#include <cmath>
#include <iostream>

using namespace std;

/// Vect2f class
///
/// 2D vector represented by its (x,y) coordinates.
/// A point really is just a vector...
class Vect2f {
public:
  Vect2f(float x = 0.0f, float y = 0.0f) : x(x), y(y){};
  float x;
  float y;
};

// Basic vector operations
Vect2f operator+(const Vect2f &a, const Vect2f &b); // a+b
Vect2f operator-(const Vect2f &a, const Vect2f &b); // a-b
Vect2f operator*(const Vect2f &a, float f);        // a*f
Vect2f operator/(const Vect2f &a, float f);        // a/f

// Dot product between two vectors
float dot(const Vect2f &a, const Vect2f &b);

// Length of a vector
float length(const Vect2f &a); // ||a||

Vect2f perpendicular(const Vect2f &a);

/// Segment class
///
/// We represent a directed line segment by its two extremities P,Q.
///
/// *--------->
/// P         Q
///
/// P,Q are 2D points represented as 2-vectors of their float
/// cartesian coordinates (x,y).
class Segment {

public:
  // Constructor
  Segment(Vect2f P_, Vect2f Q_);

  // Get coordinates of point X in the local frame defined by the Segment
  Vect2f XtoUV(Vect2f X) const;
  // Get global coordinates of point X form its coordinates u,v in local
  // frame defined by the Segment
  Vect2f UVtoX(Vect2f) const;

  // Get distance of point X from the segment
  float distance(Vect2f X) const;

  // Influence weight of this segment for spatial location X
  float weight(Vect2f X, float a, float b, float p) const;

  // Readonly accessors
  const Vect2f &getP() const { return P; };
  const Vect2f &getQ() const { return Q; };
  const Vect2f &getE1() const { return e1; };
  const Vect2f &getE2() const { return e2; };
  float getLength() const { return lPQ; };

private:
  Vect2f P;
  Vect2f Q;
  float lPQ; // length of PQ

  // We define the local frame of the segment as follows:
  //  ^
  //  |
  //
  //  e2
  //
  //  |
  //  *--- e1 --->
  //  P          Q
  //
  //  e1 = (Q-P)/(||Q-P||) in the paper
  //  e2 = perpendicalur(Q-P)/(||Q-P||) in the paper

  // Local frame
  Vect2f e1;
  Vect2f e2;
};

Image warpBy1(const Image &im, const Segment &segBefore,
              const Segment &segAfter);
Image warp(const Image &im, const vector<Segment> &segsBefore,
           const vector<Segment> &segsAfter, float a = 10.0, float b = 1.0,
           float p = 1.0);
vector<Image> morph(const Image &im1, const Image &im2,
                    const vector<Segment> &segsBefore,
                    const vector<Segment> &segsAfter, int N = 1, float a = 10.0,
                    float b = 1.0, float p = 1.0);

#endif
