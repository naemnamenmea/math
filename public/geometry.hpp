#pragma once

#include "geometry_primitives.hpp"
#include "math_constants.hpp"

enum class LINE_INTERSECTION_RESULT
{
  LINES_OVERLAP,
  LINES_ARE_PARALLEL,
  LINES_INTERSECT
};

template <class T>
double triangleArea(const point3d<T>& v1, const point3d<T>& v2) {
  return ~(v1 % v2) / 2.;
}

template <class T>
point3d<T> getPoint(const point3d<T>& beginPoint, const point3d<T>& dirVec, double segmentLen) {
  double length = ~dirVec;
  return point3d<T>(
    segmentLen * dirVec.x / length,
    segmentLen * dirVec.y / length,
    segmentLen * dirVec.z / length
    ) + beginPoint;
}

template <class T>
point3d<T> findPointInSegment(point3d<T>& m1, point3d<T>& m2, double lambda) {
  return point3d<T>(
    (m1.x + lambda * m2.x) / (1 + lambda),
    (m1.y + lambda * m2.y) / (1 + lambda),
    (m1.z + lambda * m2.z) / (1 + lambda)
    );
}

template <class T>
point3d<T> findPointOutOfSegment(point3d<T>& m0, point3d<T>& m1, double lambda) {
  return point3d<T>(
    (m0.x * (1 + lambda) - m1.x) / lambda,
    (m0.y * (1 + lambda) - m1.y) / lambda,
    (m0.z * (1 + lambda) - m1.z) / lambda
    );
}

template<class T>
double distBtw2DPointAndLine(const point3d<T>& a, const point3d<T>& b, const point3d<T>& m)
{
  line2d<T> l(a, b);
  return fabs(l.A * m.x + l.B * m.y + l.C) / sqrt(l.A * l.A + l.B * l.B);
}

template<class T>
double distBtwPointAnd2DSegment(const point3d<T>& a, const point3d<T>& b, const point3d<T>& m)
{
  if ((m - a) * (b - a) <= 0)
  {
    return ~(m - a);
  }
  else if ((m - b) * (a - b) <= 0)
  {
    return ~(m - b);
  }
  else
  {
    return distBtw2DPointAndLine(a, b, m);
  }
}

template<typename T>
double angle_between_two_vectors(const point3d<T>& a, const point3d<T>& b)
{
  return acos((a * b) / (~a) / (~b));
}

template<class T>
LINE_INTERSECTION_RESULT intersect2d(const line2d<T>& line1, const line2d<T>& line2, point3d<T>& res)
{
  double det = line1.A * line2.B - line2.A * line1.B;
  double detX = line1.C * line2.B - line2.C * line1.B;
  double detY = line1.A * line2.C - line2.A * line1.C;

  if (isZero(det))
  {
    return isZero(detX) && isZero(detY) ?
      LINE_INTERSECTION_RESULT::LINES_OVERLAP :
      LINE_INTERSECTION_RESULT::LINES_ARE_PARALLEL;
  }

  res = point3d<T>(-detX / det, -detY / det, 0);

  return LINE_INTERSECTION_RESULT::LINES_INTERSECT;
}

template<class T>
double neededPrec(T x, double __EPS) {
  T intPart = trunc(x);
  x -= intPart;
  return intPart + trunc(x / __EPS) * __EPS;
}

template<class T>
void dpcout(T x, int precision) {
  auto ss = std::cout.precision();
  std::cout << std::fixed << std::setprecision(precision) << x << std::endl;
  std::cout.precision(ss);
}