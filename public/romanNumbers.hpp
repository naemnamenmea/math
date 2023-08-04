#pragma once

#include <string_view>

namespace Converter
{
  std::string ArabicToRoman(int value);
  int RomanToArabic(std::string_view value);
}