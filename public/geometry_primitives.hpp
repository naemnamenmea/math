#pragma once

#include "point3d.hpp"

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

template <class T>
std::istream& operator>>(std::istream& is, circle<T>& c) {
  is >> c.c >> c.r;
  return is;
}

template<class T>
struct rectangle2d {
  point3d<T> a;
  point3d<T> b;
  rectangle2d() {}
  rectangle2d(const point3d<T>& a, const point3d<T>& b) : a(a), b(b) {}

  bool operator==(const rectangle2d& r) const {
    return this->a == r.a || this->b == r.b;
  }

  bool overlap2d(const rectangle2d& rect) {
    return (this->b.x > rect.a.x && rect.b.x > this->a.x || this->a.x > rect.b.x && rect.a.x > this->b.x)
      && (this->b.y > rect.a.y && rect.b.y > this->a.y || this->a.y > rect.b.y && rect.a.y > this->b.y);
  }

  bool overlap2d(const point3d<T>& l1, const point3d<T>& r1, const point3d<T>& l2, const point3d<T>& r2) {
    return !(l1.x > r2.x || l2.x > r1.x ||
      l1.y < r2.y || l2.y < r1.y);
  }
};
