#pragma once

#include <iosfwd>
#include <numeric>
#include <string>
#include <vector>

struct big_integer {
private:
  std::vector<std::uint32_t> _number;
  bool _sign = false; // false - positive, true - negative

public:
  static const std::int32_t BLOCK_SIZE = 32;
  static const std::uint32_t MAX_VALUE = std::numeric_limits<std::uint32_t>::max();
  static const std::uint64_t BASE = static_cast<std::uint64_t>(MAX_VALUE) + 1;

  std::size_t size() const noexcept;
  std::vector<std::uint32_t>& number() noexcept;
  const std::vector<std::uint32_t>& number() const noexcept;
  bool sign() const noexcept;
  std::uint32_t& operator[](std::size_t ind) noexcept;
  const std::uint32_t& operator[](std::size_t ind) const noexcept;
  void delete_zero();

public:
  big_integer();
  template <typename T>
  void constructor(T value);
  big_integer(const big_integer& other);
  big_integer(int a);
  big_integer(unsigned a);
  big_integer(long a);
  big_integer(unsigned long a);
  big_integer(long long a);
  big_integer(unsigned long long a);
  big_integer(const std::string& str);
  ~big_integer();

  void swap(big_integer& other);
  big_integer& operator=(const big_integer& other);

  std::uint32_t operator%(std::uint32_t d) const;
  big_integer& operator/=(std::int64_t d);
  std::pair<big_integer&, big_integer> divide(const big_integer& rhs);
  std::pair<big_integer&, std::uint32_t> divide(std::uint32_t d);

  big_integer& operator+=(const big_integer& rhs);
  big_integer& operator-=(const big_integer& rhs);
  big_integer& operator*=(const big_integer& rhs);
  big_integer& operator/=(const big_integer& rhs);
  big_integer& operator%=(const big_integer& rhs);

  big_integer& operator&=(const big_integer& rhs);
  big_integer& operator|=(const big_integer& rhs);
  big_integer& operator^=(const big_integer& rhs);

  big_integer& operator<<=(int rhs);
  big_integer& operator>>=(int rhs);

  big_integer operator+() const;
  big_integer operator-() const;
  big_integer operator~() const;

  big_integer& operator++();
  big_integer operator++(int);

  big_integer& operator--();
  big_integer operator--(int);

  friend big_integer operator*(const big_integer& num, std::uint32_t val);

  friend bool operator==(const big_integer& a, const big_integer& b);
  friend bool operator!=(const big_integer& a, const big_integer& b);
  friend bool operator<(const big_integer& a, const big_integer& b);
  friend bool operator>(const big_integer& a, const big_integer& b);
  friend bool operator<=(const big_integer& a, const big_integer& b);
  friend bool operator>=(const big_integer& a, const big_integer& b);

  friend std::string to_string(const big_integer& a);
};

big_integer operator+(const big_integer& a, const big_integer& b);
big_integer operator-(const big_integer& a, const big_integer& b);
big_integer operator*(const big_integer& a, const big_integer& b);
big_integer operator/(const big_integer& a, const big_integer& b);
big_integer operator%(const big_integer& a, const big_integer& b);

big_integer operator&(const big_integer& a, const big_integer& b);
big_integer operator|(const big_integer& a, const big_integer& b);
big_integer operator^(const big_integer& a, const big_integer& b);

big_integer operator<<(const big_integer& a, int b);
big_integer operator>>(const big_integer& a, int b);

bool operator==(const big_integer& a, const big_integer& b);
bool operator!=(const big_integer& a, const big_integer& b);
bool operator<(const big_integer& a, const big_integer& b);
bool operator>(const big_integer& a, const big_integer& b);
bool operator<=(const big_integer& a, const big_integer& b);
bool operator>=(const big_integer& a, const big_integer& b);

std::string to_string(const big_integer& a);
std::ostream& operator<<(std::ostream& out, const big_integer& a);
