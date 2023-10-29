
#include "morphing.h"
#include <cassert>

using namespace std;

Vect2f operator+(const Vect2f &a, const Vect2f &b) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return the vector sum of a an b
  return Vect2f(a.x+b.x, a.y+b.y); // change me
}

Vect2f operator-(const Vect2f &a, const Vect2f &b) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a-b
  return Vect2f(a.x-b.x, a.y-b.y); // change me
}

Vect2f operator*(const Vect2f &a, float f) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a*f
  return Vect2f(f*a.x, f*a.y); // change me
}

Vect2f operator/(const Vect2f &a, float f) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a/f
  return Vect2f(a.x/f, a.y/f); // change me
}

float dot(const Vect2f &a, const Vect2f &b) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return the dot product of a and b.
  return a.x*b.x+a.y*b.y; // change me
}

float length(const Vect2f &a) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return the length of a.
  return sqrtf(powf(a.x,2)+powf(a.y,2)); // change me
}

Vect2f perpendicular(const Vect2f &a) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a vector that is perpendicular to a.
  // Either direction is fine.
  return Vect2f(-a.y,a.x);
}

// The Segment constructor takes in 2 points P(x1,y1) and Q(x2,y2) corresponding
// to the ends of a segment and initialize the local reference frame e1,e2.
Segment::Segment(Vect2f P_, Vect2f Q_) : P(P_), Q(Q_) {
  // // --------- HANDOUT  PS05 ------------------------------
  // // The initializer list above ": P(P_), Q(Q_)" already copies P_
  // // and Q_, so you don't have to do it in the body of the constructor.
  // You should:
  // * Initialize the local frame e1,e2 (see header file)
  lPQ = sqrtf(powf(P.x-Q.x,2)+powf(P.y-Q.y,2));
  e1 = (Q-P)/lPQ;
  e2 = perpendicular(e1);
}

Vect2f Segment::XtoUV(Vect2f X) const {
  // --------- HANDOUT  PS05 ------------------------------
  // Compute the u,v coordinates of a point X with
  // respect to the local frame of the segment.
  // e2 ^
  //    |
  // v  +  * X
  //    | /
  //    |/
  //    *--+------>-----*
  //    P  u     e1     Q
  //                    u=1
  //
  // * Be careful with the different normalization for u and v
  return Vect2f(dot(X-P,e1), dot(X-P,e2)); // changeme
}

Vect2f Segment::UVtoX(Vect2f uv) const {
  // --------- HANDOUT  PS05 ------------------------------
  // compute the (x, y) position of a point given by the (u,v)
  // location relative to this segment.
  // * Be careful with the different normalization for u and v
  return e1*uv.x+e2*uv.y+P;
}

float Segment::distance(Vect2f X) const {
  // --------- HANDOUT  PS05 ------------------------------
  // Implement distance from a point X(x,y) to the segment. Remember the 3
  // cases from class.
  Vect2f uv = XtoUV(X);
  if (uv.x<0){
    return length(X-P);
  }else if (uv.x >lPQ){
    return length(X-Q);
  }
  return abs(uv.y);
}

Image warpBy1(const Image &im, const Segment &segBefore,
              const Segment &segAfter) {
  // --------- HANDOUT  PS05 ------------------------------
  // Warp an entire image according to a pair of segments.
 Image img = Image(im.extent(0),im.extent(1),im.extent(2));
 Vect2f t;
  for (int x=0; x<img.width();x++){
    for (int y=0;y<img.height();y++){
      t = Vect2f(x,y);
      t = segAfter.XtoUV(t);
      t = segBefore.UVtoX(t);
      for (int z=0; z<img.channels(); z++){
        img(x,y,z) = interpolateLin(im,t.x,t.y,z,true);
      }
    }
  }
  return img;
}

float Segment::weight(Vect2f X, float a, float b, float p) const {
  // --------- HANDOUT  PS05 ------------------------------
  // compute the weight of a segment to a point X(x,y) given the weight
  // parameters a,b, and p (see paper for details).
  return powf((powf(lPQ,p)/(a+distance(X))),b); // changeme
}

Image warp(const Image &im, const vector<Segment> &src_segs,
           const vector<Segment> &dst_segs, float a, float b, float p) {
  // --------- HANDOUT  PS05 ------------------------------
  // Warp an image according to a vector of before and after segments using
  // segment weighting
  Image img = Image(im.extent(0),im.extent(1),im.extent(2));
  Vect2f t,t2;
  Vect2f dsum;
  float wsum, w;
  for (int x=0; x<img.width();x++){
    for (int y=0;y<img.height();y++){
      dsum = Vect2f(0.0,0.0);
      wsum =0.0;
      t = Vect2f(x,y);
      for (int i=0; i<dst_segs.size(); i++){
        t2 = dst_segs[i].XtoUV(t);
        t2 = src_segs[i].UVtoX(t2);
        w = dst_segs[i].weight(t,a,b,p);
        dsum = dsum + (t2-t)*w;
        wsum =wsum +w;
      }
      t = t + dsum/wsum;
      for (int z=0; z<img.channels(); z++){
        img(x,y,z) = interpolateLin(im,t.x,t.y,z,true);
      }
    }
  }
  return img;
}

vector<Image> morph(const Image &im_before, const Image &im_after,
                    const vector<Segment> &segs_before,
                    const vector<Segment> &segs_after, int N, float a, float b,
                    float p) {
  // --------- HANDOUT  PS05 ------------------------------
  // return a vector of N+2 images: the two inputs plus N images that morphs
  // between im_before and im_after for the corresponding segments. im_before
  // should be the first image, im_after the last.
  if (im_before.extent(0)!=im_after.extent(0) || im_before.extent(1)!=im_after.extent(1) || im_before.extent(2)!=im_after.extent(2)){
    throw MismatchedDimensionsException();
  }
  vector<Image> imgv;
  imgv.push_back(im_before);
  Image img = Image(1,1,1);
  float t;
  for (int i=1; i<=N;i++){
    t = i/(N+1.0f);
    img = Image(im_before.extent(0),im_before.extent(1),im_before.extent(2));
    vector<Segment> tsegs;
    for (int j=0; j<segs_before.size();j++){
      Segment tseg = Segment(segs_after[j].getP()*t+segs_before[j].getP()*(1-t), segs_after[j].getQ()*t+segs_before[j].getQ()*(1-t));

      tsegs.push_back(tseg);

    }
    Image img1 = warp(im_before,segs_before,tsegs,a,b,p);
    Image img2 = warp(im_after,segs_after,tsegs,a,b,p);
    img1.write("./Output/img1morph.png");
    img2.write("./Output/img2morph.png");
    for (int x=0; x<img.width();x++){
      for (int y=0;y<img.height();y++){
        for (int z=0; z<img.channels(); z++){
          img(x,y,z) = img1(x,y,z)*(1-t) + img2(x,y,z)*t;
        }
      }
    }
    imgv.push_back(img);
  }
  imgv.push_back(im_after);
  return imgv;
}
