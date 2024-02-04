#include "big_integer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <ostream>
#include <stdexcept>

big_integer::big_integer() = default;

big_integer::big_integer(const big_integer& other) = default;

template <typename T>
void big_integer::constructor(T value) {
  if (value < 0) {
    _sign = true;
    if (value != std::numeric_limits<T>::min()) {
      value = -value;
    }
  }
  if (std::numeric_limits<T>::max() > MAX_VALUE) {
    number().push_back(value % BASE);
    number().push_back(value / BASE);
  } else {
    number().push_back(value);
  }
  delete_zero();
}

big_integer::big_integer(int a) : _sign(false) {
  constructor<int>(a);
}

big_integer::big_integer(unsigned a) : _sign(false) {
  constructor<unsigned>(a);
}

big_integer::big_integer(long long a) : _sign(false) {
  constructor<long long>(a);
}

big_integer::big_integer(unsigned long long a) : _sign(false) {
  constructor<unsigned long long>(a);
}

big_integer::big_integer(long a) : _sign(false) {
  constructor<long>(a);
}

big_integer::big_integer(unsigned long a) : _sign(false) {
  constructor<unsigned long>(a);
}

big_integer::big_integer(const std::string& str) : _sign(false) {
  if (str.empty() || str == "-") {
    throw std::invalid_argument("Invalid string");
  }
  for (std::size_t i = (str[0] == '-' ? 1 : 0); i < str.size(); i++) {
    if (!std::isdigit(static_cast<unsigned char>(str[i]))) {
      throw std::invalid_argument("Invalid string: expected digit");
    }
  }
  const std::int32_t DIGITS_COUNT = 9;
  const int STRING_BASE = 10;
  const auto MUL =
      static_cast<std::uint32_t>(std::pow(static_cast<double>(STRING_BASE), static_cast<double>(DIGITS_COUNT)));
  std::size_t last_ind = (str[0] == '-' ? 1 : 0);
  for (std::size_t i = (str[0] == '-' ? 1 : 0); i + DIGITS_COUNT < str.size(); i += DIGITS_COUNT) {
    *this *= MUL;
    *this += std::stoi(str.substr(i, DIGITS_COUNT));
    last_ind = i + DIGITS_COUNT;
  }
  if (last_ind < str.size()) {
    std::size_t len = str.size() - last_ind;
    *this *= static_cast<std::uint32_t>(std::pow(static_cast<double>(STRING_BASE), static_cast<double>(len)));
    *this += std::stoi(str.substr(last_ind, len));
  }
  if (str[0] == '-') {
    _sign = true;
  }
}

big_integer::~big_integer() = default;

std::size_t big_integer::size() const noexcept {
  return _number.size();
}

std::vector<std::uint32_t>& big_integer::number() noexcept {
  return _number;
}

const std::vector<std::uint32_t>& big_integer::number() const noexcept {
  return _number;
}

bool big_integer::sign() const noexcept {
  return _sign;
}

std::uint32_t& big_integer::operator[](std::size_t ind) noexcept {
  return number()[ind];
}

const std::uint32_t& big_integer::operator[](std::size_t ind) const noexcept {
  return number()[ind];
}

void big_integer::delete_zero() {
  while (size() > 1 && number().back() == 0) {
    number().pop_back();
  }
}

void big_integer::swap(big_integer& other) {
  std::swap(_sign, other._sign);
  std::swap(number(), other.number());
}

big_integer& big_integer::operator=(const big_integer& other) {
  if (this != &other) {
    big_integer(other).swap(*this);
  }
  return *this;
}

big_integer& big_integer::operator+=(const big_integer& rhs) {
  if (rhs.sign() != sign()) {
    _sign = !_sign;
    *this -= rhs;
    _sign = !_sign;
    return *this;
  }

  std::uint64_t carry = 0;
  number().resize(std::max(size(), rhs.size()));
  for (int i = 0; i < size(); i++) {
    std::uint64_t first = operator[](i);
    std::uint64_t second = 0;
    if (i < rhs.size()) {
      second = rhs[i];
    }

    operator[](i) = (first + second + carry) % BASE;
    carry = (first + second + carry) / BASE;
  }
  if (carry == 1) {
    number().push_back(1);
  }
  delete_zero();
  return *this;
}

bool less(const std::vector<std::uint32_t>& a, const std::vector<std::uint32_t>& b) {
  for (std::size_t i = a.size(); i > 0; i--) {
    if (a[i - 1] != b[i - 1]) {
      return a[i - 1] < b[i - 1];
    }
  }
  return false;
}

big_integer& big_integer::operator-=(const big_integer& rhs) {
  if (rhs.sign() != sign()) {
    _sign = !sign();
    *this += rhs;
    _sign = !sign();
    return *this;
  }

  bool swapped = false;
  if (size() < rhs.size() || (size() == rhs.size() && less(number(), rhs.number()))) {
    swapped = true;
    _sign = !sign();
  }

  number().resize(std::max(size(), rhs.size()));
  std::int64_t carry = 0;
  for (int i = 0; i < size(); i++) {
    std::int64_t a = operator[](i);
    std::int64_t b = 0;
    if (i < rhs.size()) {
      b = rhs[i];
    }
    if (swapped) {
      std::swap(a, b);
    }

    operator[](i) = (BASE + a - b - carry) % BASE;
    carry = ((a - b - carry >= 0) ? 0 : 1);
  }
  delete_zero();
  return *this;
}

big_integer operator*(const big_integer& num, std::uint32_t val) {
  std::uint64_t carry = 0;
  std::uint64_t a = 0;
  std::uint64_t b = val;
  big_integer result(num);
  for (std::size_t i = 0; i < num.size(); i++) {
    a = num[i];
    result[i] = (carry + a * b) % big_integer::BASE;
    carry = (carry + a * b) / big_integer::BASE;
  }
  result.number().push_back(carry);
  result.delete_zero();
  return result;
}

big_integer& big_integer::operator*=(const big_integer& rhs) {
  std::size_t this_size = size();
  number().resize(size() + rhs.size());
  for (std::size_t i = this_size; i > 0; i--) {
    std::uint64_t tmp = operator[](i - 1);
    std::uint64_t carry = 0;
    operator[](i - 1) = 0;
    for (std::size_t j = 0; j < rhs.size() || carry != 0; j++) {
      std::uint64_t current = 0;
      if (rhs.size() > j) {
        current = rhs[j];
      }
      std::uint64_t mul = current * tmp + carry;
      carry = ((operator[](i + j - 1) + mul) >> BLOCK_SIZE);
      operator[](i + j - 1) += mul;
    }
  }
  delete_zero();
  _sign = sign() != rhs.sign();
  return *this;
}

std::pair<big_integer&, std::uint32_t> big_integer::divide(std::uint32_t d) {
  std::uint64_t carry = 0;
  std::uint64_t tmp = 0;
  for (std::size_t i = size(); i > 0; i--) {
    tmp = carry * BASE + operator[](i - 1);
    operator[](i - 1) = tmp / d;
    carry = tmp % d;
  }
  delete_zero();
  return {*this, carry};
}

big_integer& big_integer::operator/=(std::int64_t d) {
  if (d < 0 || d > MAX_VALUE) {
    return *this /= big_integer(d);
  }
  return divide(d).first;
}

std::uint32_t big_integer::operator%(std::uint32_t d) const {
  std::uint64_t carry = 0;
  for (std::size_t i = size(); i > 0; i--) {
    carry = (carry * BASE + operator[](i - 1)) % d;
  }
  return carry;
}

std::pair<big_integer&, big_integer> big_integer::divide(const big_integer& rhs) {
  bool r_sign = sign();
  _sign = rhs.sign() != sign();
  if (size() < rhs.size()) {
    big_integer r = *this;
    r._sign = r_sign;
    *this = big_integer(0);
    return {*this, r};
  }
  if (rhs.size() == 1) {
    big_integer r = *this % rhs[0];
    r._sign = r_sign;
    *this /= rhs[0];
    return {*this, r};
  }

  big_integer result(0);
  big_integer r(0);
  big_integer divisor(rhs);
  divisor._sign = false;
  for (std::size_t i = size(); i > 0; i--) {
    r <<= BLOCK_SIZE;
    r += operator[](i - 1);
    result <<= BLOCK_SIZE;

    if (r >= divisor) {
      uint32_t q_guess = 0;
      uint64_t left = 0, right = BASE, mid = 0;
      while (right - left > 1) {
        mid = (left + right) / 2;
        if (divisor * mid <= r) {
          left = mid;
        } else {
          right = mid;
        }
      }
      q_guess = left;
      result += q_guess;
      r -= divisor * q_guess;
    }
  }
  result._sign = sign();
  result.delete_zero();
  r.delete_zero();
  r._sign = r_sign;
  swap(result);
  return {*this, r};
}

big_integer& big_integer::operator/=(const big_integer& rhs) {
  return *this = divide(rhs).first;
}

big_integer& big_integer::operator%=(const big_integer& rhs) {
  divide(rhs).second.swap(*this);
  return *this;
}

void not_(big_integer& bi) {
  for (auto& it : bi.number()) {
    it = ~it;
  }
  bi++;
}

big_integer& big_integer::operator&=(const big_integer& rhs) {
  bool result_sign = sign() && rhs.sign();
  number().resize(std::max(size(), rhs.size()));
  if (sign()) {
    _sign = false;
    not_(*this);
  }
  std::uint32_t carry = 0;
  std::uint32_t other = 0;
  if (rhs.sign()) {
    carry = 1;
  }
  for (std::size_t i = 0; i < size(); i++) {
    other = 0;
    if (i < rhs.size()) {
      other = rhs[i];
    }
    if (rhs.sign()) {
      other = ~other;
      if (other == MAX_VALUE && carry == 1) {
        other = 0;
        carry = 1;
      } else {
        other += carry;
        carry = 0;
      }
    }
    operator[](i) = operator[](i) & other;
  }
  delete_zero();
  if (result_sign) {
    not_(*this);
  }
  _sign = result_sign;
  delete_zero();
  return *this;
}

big_integer& big_integer::operator|=(const big_integer& rhs) {
  bool result_sign = sign() || rhs.sign();
  number().resize(std::max(size(), rhs.size()));
  if (sign()) {
    _sign = false;
    not_(*this);
  }
  std::uint32_t carry = 0;
  std::uint32_t other = 0;
  if (rhs.sign()) {
    carry = 1;
  }
  for (std::size_t i = 0; i < size(); i++) {
    other = 0;
    if (i < rhs.size()) {
      other = rhs[i];
    }
    if (rhs.sign()) {
      other = ~other;
      if (other == MAX_VALUE && carry == 1) {
        other = 0;
        carry = 1;
      } else {
        other += carry;
        carry = 0;
      }
    }
    operator[](i) = operator[](i) | other;
  }
  delete_zero();
  if (result_sign) {
    not_(*this);
  }
  _sign = result_sign;
  delete_zero();
  return *this;
}

big_integer& big_integer::operator^=(const big_integer& rhs) {
  bool result_sign = sign() ^ rhs.sign();
  number().resize(std::max(size(), rhs.size()));
  if (sign()) {
    _sign = false;
    not_(*this);
  }
  std::uint32_t carry = 0;
  std::uint32_t other = 0;
  if (rhs.sign()) {
    carry = 1;
  }
  for (std::size_t i = 0; i < size(); i++) {
    other = 0;
    if (i < rhs.size()) {
      other = rhs[i];
    }
    if (rhs.sign()) {
      other = ~other;
      if (other == MAX_VALUE && carry == 1) {
        other = 0;
        carry = 1;
      } else {
        other += carry;
        carry = 0;
      }
    }
    operator[](i) = operator[](i) ^ other;
  }
  delete_zero();
  if (result_sign) {
    not_(*this);
  }
  _sign = result_sign;
  delete_zero();
  return *this;
}

big_integer& big_integer::operator<<=(int rhs) {
  std::uint32_t zero_count = rhs / BLOCK_SIZE;
  number().resize(size() + zero_count);
  for (std::size_t i = size(); i > zero_count; i--) {
    std::swap(operator[](i - 1), operator[](i - 1 - zero_count));
  }
  std::uint64_t carry = 0;
  std::uint64_t a = 0;
  std::uint32_t shift = rhs % BLOCK_SIZE;
  for (std::size_t i = 0; i < size(); i++) {
    a = operator[](i);
    operator[](i) = (carry + (a << shift)) % BASE;
    carry = (carry + (a << shift)) / BASE;
  }
  number().push_back(carry);
  delete_zero();
  return *this;
}

big_integer& big_integer::operator>>=(int rhs) {
  std::uint32_t zero_count = rhs / BLOCK_SIZE;
  for (std::size_t i = 0; i + zero_count < size(); i++) {
    std::swap(operator[](i), operator[](i + zero_count));
  }
  for (std::size_t i = 0; i < zero_count; i++) {
    number().pop_back();
  }
  if (number().empty()) {
    number() = {0};
  }
  std::uint64_t carry = 0;
  std::uint64_t a = 0;
  std::uint32_t shift = rhs % BLOCK_SIZE;
  for (std::size_t i = size(); i > 0; i--) {
    a = carry * BASE + operator[](i - 1);
    operator[](i - 1) = (a >> shift);
    carry = a % (1 << shift);
  }
  if (sign()) {
    --(*this);
  }
  delete_zero();
  return *this;
}

big_integer big_integer::operator+() const {
  return *this;
}

big_integer big_integer::operator-() const {
  big_integer result(*this);
  result._sign = !result._sign;
  return result;
}

big_integer big_integer::operator~() const {
  return --(-*this);
}

big_integer& big_integer::operator++() {
  if (sign()) {
    _sign = !sign();
    --(*this);
    _sign = !sign();
    return *this;
  }
  bool carry = true;
  for (std::size_t i = 0; i < size(); i++) {
    if (operator[](i) == MAX_VALUE) {
      operator[](i) = 0;
    } else {
      operator[](i)++;
      carry = false;
      break;
    }
  }
  if (carry) {
    number().push_back(1);
  }
  return *this;
}

big_integer big_integer::operator++(int) {
  auto res(*this);
  operator++();
  return res;
}

big_integer& big_integer::operator--() {
  if (sign()) {
    _sign = !sign();
    ++(*this);
    _sign = !sign();
    return *this;
  }
  if (size() == 1 && operator[](0) == 0) {
    number() = {1};
    _sign = !sign();
    return *this;
  }

  for (std::size_t i = 0; i < size(); i++) {
    if (operator[](i) == 0) {
      operator[](i) = MAX_VALUE;
    } else {
      operator[](i)--;
      break;
    }
  }
  return *this;
}

big_integer big_integer::operator--(int) {
  auto res(*this);
  operator--();
  return res;
}

big_integer operator+(const big_integer& a, const big_integer& b) {
  return big_integer(a) += b;
}

big_integer operator-(const big_integer& a, const big_integer& b) {
  return big_integer(a) -= b;
}

big_integer operator*(const big_integer& a, const big_integer& b) {
  return big_integer(a) *= b;
}

big_integer operator/(const big_integer& a, const big_integer& b) {
  return big_integer(a) /= b;
}

big_integer operator%(const big_integer& a, const big_integer& b) {
  return big_integer(a) %= b;
}

big_integer operator&(const big_integer& a, const big_integer& b) {
  return big_integer(a) &= b;
}

big_integer operator|(const big_integer& a, const big_integer& b) {
  return big_integer(a) |= b;
}

big_integer operator^(const big_integer& a, const big_integer& b) {
  return big_integer(a) ^= b;
}

big_integer operator<<(const big_integer& a, int b) {
  return big_integer(a) <<= b;
}

big_integer operator>>(const big_integer& a, int b) {
  return big_integer(a) >>= b;
}

bool is_zero(const big_integer& bi) {
  return bi.size() == 0 || (bi.size() == 1 && bi[0] == 0);
}

bool operator==(const big_integer& a, const big_integer& b) {
  return (is_zero(a) && is_zero(b)) || (a.sign() == b.sign() && a.number() == b.number());
}

bool operator!=(const big_integer& a, const big_integer& b) {
  return !(a == b);
}

bool operator<(const big_integer& a, const big_integer& b) {
  if (a.size() == b.size() && a.size() == 1 && a[0] == b[0] && a[0] == 0) {
    return false;
  }
  if (a.sign() != b.sign()) {
    return a.sign() && !b.sign();
  }
  if (a.size() != b.size()) {
    return a.size() < b.size();
  }

  if (a.sign()) {
    return less(b.number(), a.number());
  } else {
    return less(a.number(), b.number());
  }
  return false;
}

bool operator>(const big_integer& a, const big_integer& b) {
  return b < a;
}

bool operator<=(const big_integer& a, const big_integer& b) {
  return !(a > b);
}

bool operator>=(const big_integer& a, const big_integer& b) {
  return !(a < b);
}

std::string to_string(const big_integer& a) {
  const uint32_t STRING_BASE = 10;
  const uint32_t DIGITS_COUNT = 9;
  const auto DIV = static_cast<uint32_t>(std::pow(static_cast<double>(STRING_BASE), static_cast<double>(DIGITS_COUNT)));
  big_integer num(a);
  std::string answer;
  while (num != big_integer(0)) {
    auto tmp = num.divide(DIV);
    std::string current = std::to_string(tmp.second);
    std::size_t new_size = answer.size() + DIGITS_COUNT;
    std::reverse(current.begin(), current.end());
    answer += current;
    if (num != big_integer()) {
      answer.resize(new_size, '0');
    }
    num.swap(tmp.first);
  }
  if (num.sign() && !answer.empty()) {
    answer.push_back('-');
  }
  std::reverse(answer.begin(), answer.end());
  if (answer.empty()) {
    answer = "0";
  }
  return answer;
}

std::ostream& operator<<(std::ostream& out, const big_integer& a) {
  return out << to_string(a);
}
