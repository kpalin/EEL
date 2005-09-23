#include <limits>
#include <math.h>


using namespace std;

#define PI 3.1415926


//returns the square of a min_k(val-k2PI,k2PI-val)
inline double squaremod(double val)
{
  double f=fabs(val);
  f-=2*PI*trunc(f/(2*PI));
  if(f>PI) {
    f=2*PI-f;
  }
#ifndef NDEBUG
  double altF=(PI-f);
  altF-=2*PI*trunc(altF/(2*PI));  // altF % (2PI)
  altF=PI-altF;


  assert(f<(altF+numeric_limits<float>::epsilon()));
  assert(f>(altF-numeric_limits<float>::epsilon()));
#endif


  assert(f<(PI+numeric_limits<float>::epsilon()));
  assert(f>(-PI-numeric_limits<float>::epsilon()));
  return f * f;
}

// The penalty score for the angle difference
inline double anglepenalty(double const d,double const D,double const nucl_per_rotation)
{
  double const theta=(d-D)*2.0*PI/nucl_per_rotation;
  return squaremod(theta)/(d+D);
}


//
// $Log$
// 
