#include "mwTPoint3d.hpp"
#include "mwMathConstants.hpp"
#include "auxGeometry.hpp"
#include "mwException.hpp"


typedef cadcam::mwTPoint3d<double> point3d;

double auxGeometry::DistBtwPointAndLine(const point3d& a, const point3d& b, const point3d& m)
{
  return 	(~((m - a) % (b - a))) / (~(b - a));
}

double auxGeometry::DistBtwPointAndSegment(const point3d& a, const point3d& b, const point3d& m)
{
  double lenAB = ~(b - a);
  double projM = (m - a) * (b - a) / lenAB;
  if (projM <= 0)
  {
    return ~(m - a);
  }
  else if (projM >= lenAB)
  {
    return ~(m - b);
  }
  else
  {
    return auxGeometry::DistBtwPointAndLine(a, b, m);
  }
}

point3d auxGeometry::PointToLineProjection(const point3d& a, const point3d& b, const point3d& m)
{
  if (cadcam::is_equal(~(a - b), 0.))
  {
    return ~(a - m) > ~(b - m) ? b : a;
  }
  point3d d(b - a);
  double t = (d * (m - a)) / (d * d);

  return  point3d(
    a.x() + d.x() * t,
    a.y() + d.y() * t,
    a.z() + d.z() * t
  );
}

point3d auxGeometry::FindClosestPointOfSegmentToPoint(const point3d& a, const point3d& b, const point3d& m)
{
  double lenAB = ~(b - a);
  double projM = (m - a) * (b - a) / lenAB;

  if (projM <= 0)
  {
    return a;
  }
  else if (projM >= lenAB)
  {
    return b;
  }
  else
  {
    return auxGeometry::PointToLineProjection(a, b, m);
  }
}

QUADRATIC_EQUATION_SOLUTION_RESULT auxGeometry::SolveQuadEq(const double a, const double b, const double c,
  double& x1, double& x2)
{
  double d = (b * b) - (4. * a * c);
  if (cadcam::is_equal(d, 0.))
  {
    x1 = -b / (2 * a);
    return QUADRATIC_EQUATION_SOLUTION_RESULT::ONE_REAL_ROOT;
  }
  else if (d > 0)
  {
    double tx1, tx2;
    tx1 = (-b + sqrt(d)) / (2 * a);
    tx2 = (-b - sqrt(d)) / (2 * a);
    x1 = (tx1 < tx2) ? tx1 : tx2;
    x2 = (tx1 > tx2) ? tx1 : tx2;
    return QUADRATIC_EQUATION_SOLUTION_RESULT::TWO_REAL_ROOTS;
  }
  else return QUADRATIC_EQUATION_SOLUTION_RESULT::IMAGINARY_ROOTS;
}

INFINITE_VERTICAL_CYLINDER_AND_LINE auxGeometry::InfiniteVerticalCylinderAndLinePosition(const point3d& cyl, const double R,
  const point3d& a, const point3d& b, std::pair<point3d, point3d>& result)
{
  auto dir = point3d(b.x() - a.x(), b.y() - a.y(), b.z() - a.z());

  double a_ = dir.x() * dir.x() + dir.y() - dir.y();
  double b_ = 2 * (dir.x() * (a.x() - cyl.x()) + dir.y() * (a.y() - cyl.y()));
  double c_ = (a.x() - cyl.x()) * (a.x() - cyl.x()) + (a.y() - cyl.y()) * (a.y() - cyl.y()) - R * R;

  double t1, t2;
  QUADRATIC_EQUATION_SOLUTION_RESULT res = auxGeometry::SolveQuadEq(a_, b_, c_, t1, t2);
  if (res == QUADRATIC_EQUATION_SOLUTION_RESULT::TWO_REAL_ROOTS)
  {
    result.first = {
      a.x() + dir.x() * t1,
      a.y() + dir.y() * t1,
      a.z() + dir.z() * t1
    };
    result.second = {
      a.x() + dir.x() * t2,
      a.y() + dir.y() * t2,
      a.z() + dir.z() * t2
    };
    return INFINITE_VERTICAL_CYLINDER_AND_LINE::INTERSECT_IN_TWO_POINTS;
  }
  else if (res == QUADRATIC_EQUATION_SOLUTION_RESULT::ONE_REAL_ROOT)
  {
    result.first = {
      a.x() + dir.x() * t1,
      a.y() + dir.y() * t1,
      a.z() + dir.z() * t1
    };
    return INFINITE_VERTICAL_CYLINDER_AND_LINE::INTERSECT_IN_ONE_POINT;
  }
  else
  {
    return INFINITE_VERTICAL_CYLINDER_AND_LINE::UNKNOWN_RESULT;
  }
}

INFINITE_VERTICAL_CYLINDER_AND_SEGMENT auxGeometry::InfiniteVerticalCylinderAndSegmentPosition(const point3d& cyl, const double R,
  const point3d& a, const point3d& b, std::pair<point3d, point3d>& result)
{
  auto res = InfiniteVerticalCylinderAndLinePosition(cyl, R, a, b, result);

  if (res == INFINITE_VERTICAL_CYLINDER_AND_LINE::INTERSECT_IN_ONE_POINT
    || res == INFINITE_VERTICAL_CYLINDER_AND_LINE::INTERSECT_IN_TWO_POINTS) {
    size_t count = 0;

    count += cadcam::is_equal(DistBtwPointAndSegment(a, b, result.first), 0.);
    count += cadcam::is_equal(DistBtwPointAndSegment(a, b, result.second), 0.);

    if (count == 2) {
      return INFINITE_VERTICAL_CYLINDER_AND_SEGMENT::INTERSECT_IN_TWO_POINTS;
    }
    else if (count == 1) {
      return INFINITE_VERTICAL_CYLINDER_AND_SEGMENT::INTERSECT_IN_ONE_POINT;
    }
    else {
      return INFINITE_VERTICAL_CYLINDER_AND_SEGMENT::UNKNOWN_RESULT;
    }
  }

  return static_cast<INFINITE_VERTICAL_CYLINDER_AND_SEGMENT>(res); // !important (the order of the values of both enumerations must match)
}



/*-----------------------
# OTHER UNUSED FEATURES #
-----------------------*/




std::pair<point3d, point3d> auxGeometry::CrossedLinesPerpendicular(
  const point3d& a_vec, const point3d& Ma,
  const point3d& b_vec, const point3d& Mb)
{
  point3d c_vec = a_vec % b_vec;

  double constXTerm = Mb.x() - Ma.x();
  double constYTerm = Mb.y() - Ma.y();
  double constZTerm = Mb.z() - Ma.z();

  double det = Determinant<double>({
    {-b_vec.x(), a_vec.x(), c_vec.x()},
    {-b_vec.y(), a_vec.y(), c_vec.y()},
    {-b_vec.z(), a_vec.z(), c_vec.z()} });

  if (cadcam::is_equal(det, 0.)) {
    throw misc::mwException(21, "lines are parallel"); // std::domain_error
  }

  double det_b_s_param = Determinant<double>({
    {constXTerm, a_vec.x(), c_vec.x()},
    {constYTerm, a_vec.y(), c_vec.y()},
    {constZTerm, a_vec.z(), c_vec.z()} });

  double det_a_t_param = Determinant<double>({
    {-b_vec.x(), constXTerm, c_vec.x()},
    {-b_vec.y(), constYTerm, c_vec.y()},
    {-b_vec.z(), constZTerm, c_vec.z()} });

  double s_param = det_b_s_param / det;
  double t_param = det_a_t_param / det;

  return {
    {a_vec.x() * t_param + Ma.x(), a_vec.y() * t_param + Ma.y(), a_vec.z() * t_param + Ma.z()},
    {b_vec.x() * s_param + Mb.x(), b_vec.y() * s_param + Mb.y(), b_vec.z() * s_param + Mb.z()}
  };
}

std::pair<point3d, point3d> auxGeometry::CrossedLinesPerpendicularQuick(const point3d& a_vec, const point3d& Ma, const point3d& b_vec, const point3d& Mb)
{
  point3d c_vec = a_vec % b_vec;

  double constXTerm = Mb.x() - Ma.x();
  double constYTerm = Mb.y() - Ma.y();
  double constZTerm = Mb.z() - Ma.z();

  double minor_1 = a_vec.y() * c_vec.z() - a_vec.z() * c_vec.y();
  double minor_2 = -b_vec.y() * c_vec.z() + b_vec.z() * c_vec.y();
  double minor_3 = -b_vec.y() * a_vec.z() + b_vec.z() * a_vec.y();
  double minor_4 = constYTerm * c_vec.z() - constZTerm * c_vec.y();
  double minor_5 = constYTerm * a_vec.z() - constZTerm * a_vec.y();
  double minor_6 = -b_vec.y() * constZTerm + b_vec.z() * constYTerm;

  double det = -b_vec.x() * minor_1 - a_vec.x() * minor_2 + c_vec.x() * minor_3;

  if (cadcam::is_equal(det, 0.)) {
    throw misc::mwException(21, "lines are parallel"); // std::domain_error
  }

  double det_b_s_param = constXTerm * minor_1 - a_vec.x() * minor_4 + c_vec.x() * minor_5;

  double det_a_t_param = -b_vec.x() * minor_4 - constXTerm * minor_2 + c_vec.x() * minor_6;

  double s_param = det_b_s_param / det;
  double t_param = det_a_t_param / det;

  return {
    {a_vec.x() * t_param + Ma.x(), a_vec.y() * t_param + Ma.y(), a_vec.z() * t_param + Ma.z()},
    {b_vec.x() * s_param + Mb.x(), b_vec.y() * s_param + Mb.y(), b_vec.z() * s_param + Mb.z()}
  };
}

bool auxGeometry::RectanglesIntersection(const point3d& aP1, const point3d& aP2, const point3d& bP1, const point3d& bP2)
{
  return (aP2.x() >= bP1.x() && bP2.x() >= aP1.x() || aP1.x() >= bP2.x() && bP1.x() >= aP2.x())
    && (aP2.y() >= bP1.y() && bP1.y() >= aP1.y() || aP1.y() >= bP1.y() && bP1.y() >= aP2.y());
}

point3d auxGeometry::PlaneAndLineIntersectionPoint(const double A, const double B, const double C, const double D,
  const point3d& s, const point3d& f)
{
  double denominator = point3d(A, B, C) * (s - f);

  if (cadcam::is_equal(point3d(A, B, C) * (s - f), 0.)) { // the line is parallel to the plane or lies in it
    throw misc::mwException(21, "divizion by zero");
  }

  double t = (point3d(A, B, C) * s + D) / denominator;

  return point3d(
    s.x() + (f.x() - s.x()) * t,
    s.y() + (f.y() - s.y()) * t,
    s.z() + (f.z() - s.z()) * t
  );
}

LINE_AND_PLANE auxGeometry::PlaneAndLinePosition(const point3d& a, const point3d& b, const point3d& c,
  const point3d& s, const point3d& f, point3d& intersectionPoint)
{
  point3d d(f - s);

  point3d n = (b - a) % (c - a);

  double A = n.x();
  double B = n.y();
  double C = n.z();
  double D = -(a * n);

  if (cadcam::is_equal(point3d(A, B, C) * (s - d) + D, 0.))
  {
    return LINE_AND_PLANE::LINE_LIES_WITHIN;
  }
  else if (cadcam::is_equal(d * n, 0.))
  {
    return LINE_AND_PLANE::PARALLEL;
  }
  else
  {
    intersectionPoint = auxGeometry::PlaneAndLineIntersectionPoint(A, B, C, D, s, f);
    return LINE_AND_PLANE::INTERSECT;
  }
}

SEGMENT_AND_PLANE auxGeometry::PlaneAndSegmentPosition(const point3d& a, const point3d& b, const point3d& c,
  const point3d& s, const point3d& f, point3d& intersectionPoint)
{
  point3d tmp;
  LINE_AND_PLANE res = auxGeometry::PlaneAndLinePosition(a, b, c, s, f, tmp);
  if (res == LINE_AND_PLANE::LINE_LIES_WITHIN)
  {
    return SEGMENT_AND_PLANE::SEGMENT_LIES_WITHIN;
  }
  else if (res == LINE_AND_PLANE::PARALLEL)
  {
    return SEGMENT_AND_PLANE::PARALLEL;
  }
  else
  {
    if (cadcam::is_equal((~(tmp - s)) + (~(f - tmp)) - (~(f - s)), 0.))
    {
      intersectionPoint = tmp;
      return SEGMENT_AND_PLANE::INTERSECT;
    }
    else
    {
      return SEGMENT_AND_PLANE::DO_NOT_INTERSECT;
    }
  }
}