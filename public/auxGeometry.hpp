#ifndef __auxGeometry_hpp__
#define __auxGeometry_hpp__

#include "mwTPoint3d.hpp"

#include <stdexcept>

enum class LINE_AND_PLANE
{
  PARALLEL,
  LINE_LIES_WITHIN,
  INTERSECT,
};

enum class SEGMENT_AND_PLANE
{
  PARALLEL,
  SEGMENT_LIES_WITHIN,
  INTERSECT,
  DO_NOT_INTERSECT,
};

enum class INFINITE_VERTICAL_CYLINDER_AND_SEGMENT
{
  INTERSECT_IN_ONE_POINT,
  INTERSECT_IN_TWO_POINTS,
  DO_NOT_INTERSECT,
  SEGMENT_LIES_WITHIN,
  UNKNOWN_RESULT,
};

enum class INFINITE_VERTICAL_CYLINDER_AND_LINE
{
  INTERSECT_IN_ONE_POINT,
  INTERSECT_IN_TWO_POINTS,
  DO_NOT_INTERSECT,
  LINE_LIES_WITHIN,
  UNKNOWN_RESULT,
};

enum class QUADRATIC_EQUATION_SOLUTION_RESULT
{
  ONE_REAL_ROOT,
  TWO_REAL_ROOTS,
  IMAGINARY_ROOTS,
};

class auxGeometry
{
  typedef cadcam::mwTPoint3d<double> point3d;
public:
  /*
  \param a - part of line
  \param b - part of line
  \param m - point
  */
  static double DistBtwPointAndLine(const point3d& a, const point3d& b, const point3d& m);
  /*
  \param a - first end of segment
  \param b - second end of segment
  \param m - point
  */
  static double DistBtwPointAndSegment(const point3d& a, const point3d& b, const point3d& m);

  static point3d PointToLineProjection(const point3d& a, const point3d& b, const point3d& m);

  static QUADRATIC_EQUATION_SOLUTION_RESULT SolveQuadEq(
    const double a, const double b, const double c,
    double& x1, double& x2);

  static point3d FindClosestPointOfSegmentToPoint(const point3d& a, const point3d& b, const point3d& m);

  static INFINITE_VERTICAL_CYLINDER_AND_LINE InfiniteVerticalCylinderAndLinePosition(
    const point3d& cyl, const double R,
    const point3d& a, const point3d& b, std::pair<point3d, point3d>& result);

  static INFINITE_VERTICAL_CYLINDER_AND_SEGMENT InfiniteVerticalCylinderAndSegmentPosition(
    const point3d& cyl, const double R,
    const point3d& a, const point3d& b, std::pair<point3d, point3d>& result);

  template<class T>
  static T Determinant(const std::vector<std::vector<T>>& Matrix)
  {
    size_t n = Matrix.size();
    if (n < 1 || Matrix[0].size() != n) {
      throw misc::mwException(21, "wrong matrix dimensions"); // std::invalid_argument
    }
    else if (n == 1)
    {
      return Matrix[0][0];
    }
    else if (n == 2)
    {
      return Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];
    }
    else
    {
      T det = 0;
      for (size_t p = 0; p < n; ++p)
      {
        std::vector<std::vector<T>> CofactorMatrix(n - 1, std::vector<T>(n - 1));
        for (size_t i = 1; i < n; ++i)
        {
          size_t cofactor_j = 0;
          for (size_t j = 0; j < n; ++j)
          {
            if (j == p) continue;

            CofactorMatrix[i - 1][cofactor_j] = Matrix[i][j];
            cofactor_j++;
          }
        }
        det += static_cast<T>(pow(-1, p)) * Matrix[0][p] * Determinant(move(CofactorMatrix));
      }
      return det;
    }
  }

  static std::pair<point3d, point3d> CrossedLinesPerpendicular(
    const point3d& a_vec, const point3d& Ma,
    const point3d& b_vec, const point3d& Mb);

  static std::pair<point3d, point3d> CrossedLinesPerpendicularQuick(
    const point3d& a_vec, const point3d& Ma,
    const point3d& b_vec, const point3d& Mb);




  /*-----------------------
  # OTHER UNUSED FEATURES #
  -----------------------*/




  static bool RectanglesIntersection(const point3d& aP1, const point3d& aP2, const point3d& bP1, const point3d& bP2);
  /*
  \param coefficients A,B,C,D determines plane
  \param 3dpoints s,f determines line
  */
  static point3d PlaneAndLineIntersectionPoint(
    const double A, const double B, const double C, const double D,
    const point3d& s, const point3d& f);
  /*
  \param points a,b,c determines plane
  \param 3dpoints s,f determines line
  \output intersectionPoint if it exist
  */
  static LINE_AND_PLANE PlaneAndLinePosition(
    const point3d& a, const point3d& b, const point3d& c,
    const point3d& s, const point3d& f, point3d& intersectionPoint);
  /*
  \param points a,b,c determines plane
  \param 3dpoints s,f determines segment
  \output intersectionPoint if it exist
  */
  static SEGMENT_AND_PLANE PlaneAndSegmentPosition(
    const point3d& a, const point3d& b, const point3d& c,
    const point3d& s, const point3d& f, point3d& intersectionPoint);
protected:
private:
};

#endif //__auxGeometry_hpp__