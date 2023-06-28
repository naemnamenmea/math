#define _USE_MATH_DEFINES

#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <array>
#include <map>
#include <iomanip>
#include <cmath>

#ifdef LOCAL_LAUNCH
#include <fstream>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;
#endif

constexpr auto __EPS = 1e-9;
constexpr auto __MOD = 1e9 + 7;

enum class LINE_INTERSECTION_RESULT
{
  LINES_OVERLAP,
  LINES_ARE_PARALLEL,
  LINES_INTERSECT
};

template<class T>
T gcd(T a, T b) {
  return b == 0 ? a : gcd(b, a % b);
}

template <typename T>
inline bool isZero(T val)
{
  return (val > 0 ? val : -val) < __EPS;
}

template <typename T>
inline bool isEqual(T a, T b)
{
  return isZero(a - b);
}

template<class T>
struct point3d;

template <class T = double>
struct point2d {
  double x;
  double y;

  point2d<T>& operator=(const point3d<T>& p) {
    x = p.x;
    y = p.y;

    return *this;
  }

  inline bool operator==(const point2d<T>& p) const
  {
    return this->x == p.x && this->y == p.y;
  }
};

template <class T = double>
struct point3d {

  point3d(T x = 0, T y = 0, T z = 0) : x(x), y(y), z(z) {}

  point3d(const point2d<T>& p) : x(p.x), y(p.y), z(0) {}

  point3d& operator=(const point2d<T>& p) {
    x = p.x;
    y = p.y;
    z = 0;

    return *this;
  }

  inline T operator*(const point3d<T>& p) const
  {
    return this->x * p.x + this->y * p.y + this->z * p.z;
  }

  inline point3d<T> operator-(const point3d<T>& p) const
  {
    return point3d<T>(this->x - p.x, this->y - p.y, this->z - p.z);
  }

  inline point3d<T> operator+(const point3d<T>& p) const
  {
    return point3d<T>(this->x + p.x, this->y + p.y, this->z + p.z);
  }

  inline point3d<T> operator/(const double factor) const
  {
    return point3d<T>(this->x / factor, this->y / factor, this->z / factor);
  }

  inline point3d<T> operator*(const double factor) const
  {
    return point3d<T>(this->x * factor, this->y * factor, this->z * factor);
  }

  inline bool operator==(const point3d<T>& p) const
  {
    return this->x == p.x && this->y == p.y && this->z == p.z;
  }

  inline point3d<T> operator-() const
  {
    return point3d<T>(-this->x, -this->y, -this->z);
  }

  inline void operator+=(const point3d<T>& p) {
    this->x += p.x;
    this->y += p.y;
    this->z += p.z;
  }

  inline double operator~() const
  {
    return sqrt((*this) * (*this));
  }

  point3d<T> operator-(const T v) const {
    return point3d<T>(this->x - v, this->y - v, this->z - v);
  }
  point3d<T> operator+(const T v) const {
    return point3d<T>(this->x + v, this->y + v, this->z + v);
  }

  bool operator!=(const point3d<T>& p1) const {
    return !(*this == p1);
  }

  point3d<T> operator%(const point3d<T>& p) const {
    return point3d<T>(
      this->y * p.z - this->z * p.y,
      -(this->x * p.z - this->z * p.x),
      this->x * p.y - this->y * p.x
      );
  }

  T x;
  T y;
  T z;
};

template<class T>
struct surface {
  surface() {}
  surface(T A, T B, T C, T D) : A(A), B(B), C(C), D(D) {}

  T A, B, C, D;
};

template<class T>
struct sphere {
  sphere() {}
  sphere(const point3d<T>& c, double r) : c(c), r(r) {}

  bool intersect(const sphere& s) const {
    return !(~(this->c - s.c) >= this->r + s.r);
  }

  point3d<T> c;
  double r;
};

template <class T>
struct line2d {
  line2d() {}
  line2d(const point3d<T>& a, const point3d<T>& b) : A(a.y - b.y), B(-(a.x - b.x)), C(a.x* b.y - b.x * a.y) {}
  // line2d(const point3d<T>& dot, const point3d<T>& normal) :A(normal.x), B(normal.y), C(-dot * normal) {}
  line2d(double a, double b, double c) : A(a), B(b), C(c) {}

  bool evY(double x, double y) {
    return isZero(A * x + B * y + C);
  }

  T A;
  T B;
  T C;
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

template <class T>
struct circle {
  circle() {}
  circle(const point3d<T>& c, double r) :c(c), r(r) {}

  double area() const {
    return M_PI * this->r * this->r;
  }

  double dualArea(const circle& c) const {
    double d = ~(this->c - c.c);
    if (d >= this->r + c.r) {
      return this->area() + c.area();
    }
    else if (d <= fabs(this->r - c.r)) {
      return (this->r > c.r ? this->area() : c.area());
    }
    else {
      double f1 = 2 * acos((this->r * this->r - c.r * c.r + d * d) / (2 * this->r * d));
      double f2 = 2 * acos((c.r * c.r - this->r * this->r + d * d) / (2 * c.r * d));
      double s1 = this->r * this->r * (f1 - sin(f1)) / 2.;
      double s2 = c.r * c.r * (f2 - sin(f2)) / 2.;
      //debug(s1, s2);
      return (this->area() + c.area()) - (s1 + s2);
    }
  }

  point3d<T> c;
  double r;
};

template<class T>
struct rectangle2d {
  point3d<T> a;
  point3d<T> b;
  rectangle2d() {}
  rectangle2d(const point3d<T>& a, const point3d<T>& b) : a(a), b(b) {}

  bool operator==(const rectangle2d<T>& r) const {
    return this->a == r.a || this->b == r.b;
  }

  bool overlap2d(const rectangle2d<T>& rect) {
    return (this->b.x > rect.a.x && rect.b.x > this->a.x || this->a.x > rect.b.x && rect.a.x > this->b.x)
      && (this->b.y > rect.a.y && rect.b.y > this->a.y || this->a.y > rect.b.y && rect.a.y > this->b.y);
  }

  bool overlap2d(const point3d<T>& l1, const point3d<T>& r1, const point3d<T>& l2, const point3d<T>& r2) {
    return !(l1.x > r2.x || l2.x > r1.x ||
      l1.y < r2.y || l2.y < r1.y);
  }

};

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

template <class T>
std::ostream& operator<<(std::ostream& os, const point3d<T>& p) {
  os << p.x << ' ' << p.y << ' ' << p.z;
  return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const point2d<T>& p) {
  os << p.x << ' ' << p.y;
  return os;
}

template <class T>
std::istream& operator>>(std::istream& is, point3d<T>& p) {
  is >> p.x >> p.y >> p.z;
  return is;
}

template <class T>
std::istream& operator>>(std::istream& is, point2d<T>& p) {
  is >> p.x >> p.y;
  return is;
}

template <class T>
std::istream& operator>>(std::istream& is, circle<T>& c) {
  is >> c.c >> c.r;
  return is;
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

using namespace std;

template<class T>
void normalize(T& a, T& b) {
  T gcd_v = gcd(a, b);
  a /= gcd_v;
  b /= gcd_v;

  if (b < 0) {
    b = -b;
    a = -a;
  }
}

int main()
{
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

#ifdef LOCAL_LAUNCH
  fs::path root = fs::path(__FILE__).parent_path();
  ifstream ifs(root / "input.txt"); cin.rdbuf(ifs.rdbuf());
#endif

  

  return 0;
}