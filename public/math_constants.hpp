#pragma once

#include <cmath>
#include <limits>
#ifdef _MSC_VER
#include <immintrin.h>
#endif

inline unsigned int MostSignificantBitIdx(unsigned long long x)
{
  if (x == 0) return static_cast<unsigned int>(-1);
#ifdef _MSC_VER
  return static_cast<unsigned int>(8 * sizeof(unsigned long long) - _lzcnt_u64(x) - 1);
#else
  return static_cast<unsigned int>(8 * sizeof(unsigned long long) - __builtin_clzll(x) - 1);
#endif
}

namespace mathdef
{
  const double MATH_TOL = 1e-12;
  const float MATH_TOL_FLOAT = 1e-6f;
  const double MATH_TOL_LD = 1e-16;

  const double PI = 3.1415926535;
  const double D2R = PI / 180;

  constexpr double MAX_DOUBLE = std::numeric_limits<double>::max();

  inline bool is_within(
    const size_t value,
    const size_t from,
    const size_t to,
    const bool isIncludeBegin = true,
    const bool isIncludeEnd = true
  )
  {
    if (!(isIncludeBegin ? from <= value : from < value))
      return false;
    if (!(isIncludeEnd ? value <= to : value < to))
      return false;
    return true;
  }

  inline double math_tol(const double&)
  {
    return MATH_TOL;
  }

  inline bool is_eq(const double& a, const double& b)
  {
    return std::abs(a - b) < MATH_TOL;
  }

  inline bool is_neq(const double& val1, const double& val2)
  {
    return !is_eq(val1, val2);
  }

  inline bool is_eq(const long double& a, const long double& b)
  {
    return std::abs(a - b) < MATH_TOL_LD;
  }

  inline bool is_eq(const double& a, const double& b, const double& tol)
  {
    return std::abs(a - b) < tol;
  }

  inline bool is_lt(const double& toCompare, const double& source)
  {
    return toCompare < (source - MATH_TOL);
  }

  inline bool is_gt(const double& toCompare, const double& source)
  {
    return toCompare > (source + MATH_TOL);
  }

  inline bool is_lte(const double& toCompare, const double& source)
  {
    return toCompare <= (source + MATH_TOL);
  }

  inline bool is_gte(const double& toCompare, const double& source)
  {
    return toCompare >= (source - MATH_TOL);
  }

  inline bool is_lte(const long double& toCompare, const long double& source)
  {
    return toCompare <= (source + MATH_TOL);
  }

  inline bool is_zero(const long double& value)
  {
    return std::abs(value) <= MATH_TOL_LD;
  }

  inline bool is_not_zero(const long double& value)
  {
    return !is_zero(value);
  }
}  // namespace mathdef
