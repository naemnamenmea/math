#include "romanNumbers.hpp"
#include <unordered_map>
#include <stdexcept>
#include <cmath>

using namespace std;

namespace
{
  const unordered_map<char, int> m{ {'I',1},{'V',5},{'X',10},{'L',50},{'C',100},{'D',500},{'M',1000} };

  int GetLength(int num)
  {
    int res = 0;
    while (num > 0)
    {
      ++res;
      num /= 10;
    }

    return res;
  }

  char GetPowerChar1(int power)
  {
    char res;
    switch (power)
    {
    case 3: res = 'M'; break;
    case 2: res = 'C'; break;
    case 1: res = 'X'; break;
    case 0: res = 'I'; break;
    default: throw std::runtime_error("");
    }

    return res;
  }

  char GetPowerChar5(int power)
  {
    char res;
    switch (power)
    {
    case 2: res = 'D'; break;
    case 1: res = 'L'; break;
    case 0: res = 'V'; break;
    default: throw std::runtime_error("");
    }

    return res;
  }
}

namespace Converter
{
  string ArabicToRoman(int num)
  {
    if (num < 1)
      throw runtime_error("number is out of roman number bounds");

    string res;

    int len = GetLength(num);
    for (int i = 0; i < len; ++i)
    {
      int power = len - i - 1;
      int digit = (num / (int)std::pow(10, power)) % 10;

      char c = GetPowerChar1(power);
      if (digit < 4)
      {
        res += std::string(digit, c);
      }
      else if (digit == 4)
      {
        char c2 = GetPowerChar5(power);
        res += c;
        res += c2;
      }
      else if (digit < 9)
      {
        char c2 = GetPowerChar5(power);
        res += c2 + std::string(digit - 5, c);
      }
      else
      {
        char c2 = GetPowerChar1(power + 1);
        res += c;
        res += c2;
      }
    }

    return res;
  }

  int RomanToArabic(string_view s)
  {
    if (s.empty())
      throw runtime_error("Invalid roman number");

    int res = 0;
    int prev = (int)1e6;

    for (char c : s)
    {
      int cur = m.at(c);
      if (cur > prev) {
        res -= 2 * prev;
      }
      res += cur;
      prev = cur;
    }

    return res;
  }
}