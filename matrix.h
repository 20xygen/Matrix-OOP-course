#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <type_traits>
#include <cassert>
#include <array>
#include <ctime>

#define CPP23 1


// 'biginteger.h' {

enum class Sign : int {
  Negative = -1,
  Neutral = 0,
  Positive = 1
};

Sign operator-(const Sign sign) {
  return static_cast<Sign>(-static_cast<int>(sign));
}

class BigInteger {
private:
  static const int kTo_string_base_ = 10;
  static const int kBase_ = 1e9;
  static const int kExp_ = 9;
  std::vector<int> bits_;
  Sign sign_;

  static void CutZeros(std::vector<int>& input) {
    if (input.back() != 0) {
      return;
    }
    int pin = input.size() - 1;
    while (pin > 0 && input[pin] == 0) {
      input.pop_back();
      --pin;
    }
  }

  void Fit() {
    for (int pos = bits_.size() - 1; pos > 0; --pos) {
      if (bits_[pos] == 0) {
        bits_.pop_back();
      } else {
        break;
      }
    }
    if (bits_.size() == 1 && bits_[0] == 0) {
      sign_ = Sign::Neutral;
    }
  }

  static void Add(std::vector<int>& first, const std::vector<int>& second) {
    std::vector<int> neutral = {0};
    if (first == neutral) {
      first = second;
      return;
    }
    if (second == neutral) {
      return;
    }
    while (first.size() < second.size()) {
      first.push_back(0);
    }
    first.push_back(0);
    int buf = 0;
    for (size_t i = 0; i < second.size(); ++i) {
      first[i] += second[i] + buf;
      if (first[i] >= kBase_) {
        first[i] -= kBase_;
        buf = 1;
      } else {
        buf = 0;
      }
    }
    if (buf > 0) {
      first[second.size()] += buf;
    }
  }

  static void Add(const std::vector<int>& first, const std::vector<int>& second, std::vector<int>& destination) {
    destination = first;
    Add(destination, second);
  }

  static void Subtract(std::vector<int>& first, const std::vector<int>& second) {
    const std::vector<int> third = first;
    Subtract(third, second, first);
  }

  static void Subtract(const std::vector<int>& first, const std::vector<int>& second, std::vector<int>& destination) {
    std::vector<int> neutral = {0};
    if (first == neutral) {
      destination = second;
      return;
    }
    if (second == neutral) {
      destination = first;
      return;
    }
    destination = {};
    while (destination.size() < first.size()) {
      destination.push_back(0);
    }
    int buf = 0;
    for (size_t i = 0; i < second.size(); ++i) {
      if (first[i] - buf < second[i]) {
        destination[i] = first[i] + kBase_ - buf - second[i];
        buf = 1;
      } else {
        destination[i] = first[i] - buf - second[i];
        buf = 0;
      }
    }
    for (size_t i = second.size(); i < first.size(); ++i) {
      if (first[i] - buf < 0) {
        destination[i] = first[i] + kBase_ - buf;
        buf = 1;
      } else {
        destination[i] = first[i] - buf;
        buf = 0;
      }
    }
  }

  static void InversedSubtract(std::vector<int>& first, const std::vector<int>& second) {
    std::vector<int> neutral = {0};
    if (first == neutral) {
      first = second;
      return;
    }
    if (second == neutral) {
      return;
    }
    while (first.size() < second.size()) {
      first.push_back(0);
    }
    int buf = 0;
    for (size_t i = 0; i < second.size(); ++i) {
      if (second[i] - buf < first[i]) {
        first[i] = second[i] + kBase_ - buf - first[i];
        buf = 1;
      } else {
        first[i] = second[i] - buf - first[i];
        buf = 0;
      }
    }
    if (buf > 0) {
      first[second.size()] -= buf;
    }
  }

  static void Multiply(const std::vector<int>& first, const std::vector<int>& second, std::vector<int>& destination) {
    std::vector<int> neutral = {1};
    if (first == neutral) {
      destination = second;
      return;
    }
    if (second == neutral) {
      destination = first;
      return;
    }
    std::vector<int> deleter = {0};
    if (first == deleter || second == deleter) {
      destination = deleter;
      return;
    }
    destination = std::vector<int>(first.size() + second.size() + 1, 0);
    std::vector<long long> destiny = std::vector<long long>(first.size() + second.size() + 1, 0);
    for (size_t i = 0; i < second.size(); ++i) {
      long long buf = 0;
      for (size_t j = 0; j < first.size(); ++j) {
        long long nu = 1ll * first[j] * second[i] + buf;
        buf = nu / kBase_;
        destiny[i + j] += nu % kBase_;
      }
      destiny[i + first.size()] += buf;
    }
    long long buf = 0;
    for (size_t i = 0; i < destiny.size(); ++i) {
      int bu = buf;
      buf = (destiny[i] + buf) / kBase_;
      destiny[i] = (destiny[i] + bu) % kBase_;
    }
    for (size_t i = 0; i < destiny.size(); ++i) {
      destination[i] = destiny[i];
    }
  }

  static bool NotHigher(const std::vector<int>& first, const std::vector<int>& second) {
    int delta_1 = 0;
    int delta_2 = 0;
    while (first[first.size() - 1 - delta_1] == 0 && first.size() - 1 - delta_1 > 0) {
      ++delta_1;
    }
    while (second[second.size() - 1 - delta_2] == 0 && second.size() - 1 - delta_2 > 0) {
      ++delta_2;
    }
    if (first.size() - delta_1 < second.size() - delta_2) {
      return true;
    }
    if (first.size() - delta_1 > second.size() - delta_2) {
      return false;
    }
    for (ssize_t i = first.size() - delta_1 - 1; i >= 0; --i) {
      if (first[i] > second[i]) {
        return false;
      }
      if (first[i] < second[i]) {
        return true;
      }
    }
    return true;
  }

  static void Divide(const std::vector<int>& first, const std::vector<int>& second, std::vector<int>& destination) {
    std::vector<int> neutral = {1};
    if (second == neutral) {
      destination = first;
      return;
    }
    std::vector<int> deleter = {0};
    if (first == deleter) {
      destination = deleter;
      return;
    }
    destination = std::vector<int>(first.size() + 1, -1);
    int pin = first.size();
    std::vector<int> buf;
    for (int i = first.size() - 1; i >= 0; --i) {
      buf.insert(buf.begin(), first[i]);
      for (int pos = buf.size() - 1; pos > 0; --pos) {
        if (buf[pos] == 0) {
          buf.pop_back();
        } else {
          break;
        }
      }
      std::vector<int> buf_2 = {0};
      int left = 0;
      int right = kBase_ - 1;
      int mid;
      while (right > left + 1) {
        mid = left + (right - left) / 2;
        buf_2 = {};
        std::vector<int> mi;
        mi.push_back(mid);
        Multiply(second, mi, buf_2);
        if (NotHigher(buf, buf_2)) {
          right = mid;
        } else {
          left = mid;
        }
      }
      buf_2 = {};
      Multiply(std::vector<int>({right}), second, buf_2);
      CutZeros(buf_2);
      int number;
      if (NotHigher(buf_2, buf)) {
        number = right;
      } else {
        Subtract(buf_2, second);
        CutZeros(buf_2);
        number = left;
      }
      destination[pin] = number;
      --pin;
      Subtract(buf, buf_2);
    }
    int start = pin + 1;
    for (size_t i = start; i < destination.size(); ++i) {
      destination[i - start] = destination[i];
    }
    for (size_t i = destination.size() - start; i < destination.size(); ++i) {
      destination.pop_back();
    }
  }

  BigInteger(std::vector<int>& source) {
    bits_ = source;
    sign_ = Sign::Positive;
  }

public:
  BigInteger(int integer) {
    if (integer > 0) {
      sign_ = Sign::Positive;
    } else if (integer < 0) {
      sign_ = Sign::Negative;
    } else {
      sign_ = Sign::Neutral;
    }
    if (integer < 0) {
      integer *= -1;
    }
    if (integer < kBase_) {
      bits_.push_back(integer);
      return;
    }
    if (integer == kBase_) {
      bits_.push_back(1);
      bits_.push_back(0);
      return;
    }
    int max = 0;
    long divider = 1;
    while (divider * kBase_ <= integer) {
      ++max;
      divider *= kBase_;
    }
    divider = kBase_;
    for (; max >= 0; --max) {
      bits_.push_back((integer % divider) / (divider / kBase_));
      divider *= kBase_;
    }
  }

  BigInteger() {
    sign_ = Sign::Neutral;
    bits_.push_back(0);
  };

  BigInteger(std::string str) {
    if (str[0] == '-') {
      sign_ = Sign::Negative;
      str.erase(str.begin());
    } else if (str == "0") {
      sign_ = Sign::Neutral;
      bits_ = {0};
      return;
    } else {
      sign_ = Sign::Positive;
    }
    std::string buf;
    ssize_t i;
    for (i = str.size() - 1; i >= kExp_ - 1; i -= kExp_) {
      buf = "";
      for (int j = kExp_ - 1; j >= 0; j--) {
        buf += str[i - j];
      }
      bits_.push_back(std::atoi(buf.c_str()));
    }
    buf = "";
    for (int j = 0; j <= i; ++j) {
      buf += str[j];
    }
    if (buf.size() > 0) {
      bits_.push_back(std::atoi(buf.c_str()));
    }
  }

  BigInteger(const char* str) {
    std::string st = str;
    *this = BigInteger(st);
  }

  std::string toString() const {
    std::string str;
    for (size_t i = 0; i < bits_.size() ; ++i) {
      int temp = bits_[i];
      for (int j = 0; j < kExp_; ++ j) {
        char ch = temp % kTo_string_base_ + '0';
        temp /= kTo_string_base_;
        str.push_back(ch);
      }
    }
    while (str.back() == '0' && str.size() > 1) {
      str.pop_back();
    }
    if (sign_ == Sign::Negative) {
      str.push_back('-');
    }
    str = {str.rbegin(), str.rend()};
    return str;
  }

  BigInteger operator-() const {
    BigInteger bigInteger = BigInteger(*this);
    bigInteger.sign_ = -sign_;
    return bigInteger;
  }

  friend bool operator==(const BigInteger&, const BigInteger&);

  friend bool operator<(const BigInteger&, const BigInteger&);

  BigInteger& operator+=(const BigInteger& second) {
    if (sign_ == Sign::Neutral) {
      *this = second;
      return *this;
    }
    if (second.sign_ == Sign::Neutral) {
      return *this;
    }
    if (sign_ == second.sign_) {
      Add(bits_, second.bits_);
      Fit();
      return *this;
    }
    if (NotHigher(bits_, second.bits_)) {
      InversedSubtract(bits_, second.bits_);
      sign_ = second.sign_;
      Fit();
      return *this;
    }
    Subtract(bits_, second.bits_);
    Fit();
    return *this;
  }

  BigInteger& operator-=(const BigInteger& second) {
    if (sign_ == Sign::Neutral) {
      *this = -second;
      return *this;
    }
    if (second.sign_ == Sign::Neutral) {
      return *this;
    }
    if (sign_ != second.sign_) {
      Add(bits_, second.bits_);
      Fit();
      return *this;
    }
    if (*this == second) {
      *this = 0;
      return *this;
    }
    if (NotHigher(bits_, second.bits_)) {
      InversedSubtract(bits_, second.bits_);
      sign_ = -second.sign_;
      Fit();
      return *this;
    }
    Subtract(bits_, second.bits_);
    Fit();
    return *this;
  }

  BigInteger& operator*=(const BigInteger& second) {
    if (sign_ == Sign::Neutral) {
      return *this;
    }
    if (second.sign_ == Sign::Neutral) {
      *this = second;
      return *this;
    }
    BigInteger created;
    std::vector<int> source;
    Multiply(bits_, second.bits_, source);
    created = BigInteger(source);
    if (sign_ != second.sign_) {
      created.sign_ = Sign::Negative;
    }
    *this = created;
    Fit();
    return *this;
  }

  BigInteger& operator/=(const BigInteger& second) {
    assert(second.sign_ != Sign::Neutral);
    if (*this == BigInteger(1)) {
      *this = second;
      return *this;
    }
    if (*this == -BigInteger(1)) {
      *this = -second;
      return *this;
    }
    if (*this == BigInteger(0)) {
      *this = BigInteger(0);
      return *this;
    }
    if (second == *this) {
      *this = BigInteger(1);
      return *this;
    }
    if (second == -(*this)) {
      *this = -BigInteger(1);
      return *this;
    }
    if (NotHigher(bits_, second.bits_)) {
      *this = BigInteger(0);
      return *this;
    }
    std::vector<int> source;
    Divide(bits_, second.bits_, source);
    if (sign_ != second.sign_) {
      *this = BigInteger(source);
      sign_ = Sign::Negative;
    } else {
      *this = BigInteger(source);
    }
    Fit();
    return *this;
  }

  BigInteger& operator++() {
    *this += BigInteger(1);
    return *this;
  }

  BigInteger& operator--() {
    *this -= BigInteger(1);
    return *this;
  }

  BigInteger operator++(int) {
    BigInteger bigInteger = *this;
    *this += BigInteger(1);
    return bigInteger;
  }

  BigInteger operator--(int) {
    BigInteger bigInteger = *this;
    *this -= BigInteger(1);
    return bigInteger;
  }

  explicit operator bool() const { return sign_ != Sign::Neutral; }

  explicit operator int() const {
    return atoi(toString().c_str());
  }

  BigInteger& operator%=(const BigInteger&);
};

std::ostream& operator<<(std::ostream& out, const BigInteger& biginteger) {
  return out << biginteger.toString();
}

std::istream& operator>>(std::istream& in, BigInteger& biginteger) {
  std::string str;
  in >> str;
  biginteger = BigInteger(str);
  return in;
}

bool operator==(const BigInteger& first, const BigInteger& second) {
  if (second.bits_.size() != first.bits_.size() || second.sign_ != first.sign_) {
    return false;
  }
  for (size_t i = 0; i < first.bits_.size(); ++i) {
    if (second.bits_[i] != first.bits_[i]) {
      return false;
    }
  }
  return true;
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
  return !(first == second);
}

bool operator<(const BigInteger& first, const BigInteger& second) {
  if (first.sign_ < second.sign_) {
    return true;
  }
  if (first.sign_ > second.sign_) {
    return false;
  }
  if (first.sign_ == Sign::Neutral) {
    return false;
  }
  bool modifier = first.sign_ > Sign::Neutral;
  if (first.bits_.size() < second.bits_.size()) {
    return modifier;
  }
  if (first.bits_.size() > second.bits_.size()) {
    return !modifier;
  }
  for (ssize_t i = first.bits_.size() - 1; i >= 0; --i) {
    if (first.bits_[i] > second.bits_[i]) {
      return !modifier;
    }
    if (first.bits_[i] < second.bits_[i]) {
      return modifier;
    }
  }
  return false;
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
  return !(first < second);
}

bool operator>(const BigInteger& first, const BigInteger& second) {
  return second < first;
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
  return !(first > second);
}

BigInteger operator"" _bi(const char* number) {
  return BigInteger(number);
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
  BigInteger ret = first;
  ret += second;
  return ret;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
  BigInteger ret = first;
  ret -= second;
  return ret;
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
  BigInteger ret = first;
  ret *= second;
  return ret;
}

BigInteger operator/(const BigInteger& first, const BigInteger& second) {
  BigInteger ret = first;
  ret /= second;
  return ret;
}

BigInteger& BigInteger::operator%=(const BigInteger& second) {
  assert(second.sign_ != Sign::Neutral);
  *this =  *this - (*this / second) * second;
  Fit();
  return *this;
}

BigInteger operator%(const BigInteger& first, const BigInteger& second) {
  BigInteger ret = first;
  ret %= second;
  return ret;
}

BigInteger GCD(BigInteger first, BigInteger second) {
  if (first < 0) {
    first = -first;
  }
  if (second < 0) {
    second = -second;
  }
  while(first != 0 && second != 0) {
    if (first > second) {
      first %= second;
      continue;
    }
    second %= first;
  }
  return (first == 0) ? second : first;
}

class Rational {
private:
  static const int kDouble_precision_ = 20;
  BigInteger numerator_;
  BigInteger denominator_;

  Rational(const BigInteger& numerator, const BigInteger& denominator){
    BigInteger zero = 0;
    if (denominator < zero) {
      denominator_ = -denominator;
      numerator_ = -numerator;
    } else {
      denominator_ = denominator;
      numerator_ = numerator;
    }
  }

  void Fit() {
    BigInteger core = GCD(numerator_, denominator_);
    numerator_ /= core;
    denominator_ /= core;
  }

public:
  Rational() {
    numerator_ = 0;
    denominator_ = 1;
  }

  Rational(const BigInteger& bigInteger) {
    numerator_ = bigInteger;
    denominator_ = 1;
  }

  Rational(int integer) {
    numerator_ = integer;
    denominator_ = 1;
  }

  Rational operator-() const {
    Rational created = Rational(-numerator_, denominator_);
    return created;
  }

  friend bool operator==(const Rational&, const Rational&);

  friend bool operator<(const Rational&, const Rational&);

  Rational& operator-=(const Rational& second) {
    BigInteger core = GCD(denominator_, second.denominator_);
    *this = Rational(numerator_ * (second.denominator_ / core) - second.numerator_ * (denominator_ / core), denominator_ * second.denominator_ / core);
    Fit();
    return *this;
  }

  Rational& operator+=(const Rational& second) {
    *this -= -second;
    return *this;
  }

  Rational& operator*=(const Rational& second) {
    *this = Rational(numerator_ * second.numerator_, denominator_ * second.denominator_);
    Fit();
    return *this;
  }

  Rational& operator/=(const Rational& second) {
    assert(second.numerator_ != 0);
    *this = Rational(numerator_ * second.denominator_, denominator_ * second.numerator_);
    Fit();
    return *this;
  }

  std::string toString() const {
    std::string str = numerator_.toString();
    if (denominator_ != BigInteger(1)) {
      str += "/" + denominator_.toString();
    }
    return str;
  }

  std::string asDecimal(size_t precision = 0) const {
    assert(precision >= 0);
    std::string str = (numerator_ / denominator_).toString();
    if (str == "0" && numerator_ < 0 && precision > 0) {
      str = "-0";
    }
    if (precision == 0) {
      return str;
    }
    str += ".";
    BigInteger num;
    if (numerator_ < 0) {
      num = -numerator_;
    } else {
      num = numerator_;
    }
    num %= denominator_;
    BigInteger factor = 1;
    for (size_t i = 0; i < precision; ++i) {
      factor *= 10;
    }
    num *= factor;
    std::string small = (num / denominator_).toString();
    if (precision - small.size() > 0) {
      str += std::string(precision - small.size(), '0');
    }
    str += small;
    return str;
  }

  explicit operator double() const {
    std::string str = asDecimal(kDouble_precision_);
    return std::stod(str);
  }
};

bool operator==(const Rational& first, const Rational& second) {
  return first.numerator_ == second.numerator_ && first.denominator_ == second.denominator_;
}

bool operator!=(const Rational& first, const Rational& second) {
  return !(first == second);
}

bool operator<(const Rational& first, const Rational& second) {
  BigInteger core = GCD(first.denominator_, second.denominator_);
  return first.numerator_ * (second.denominator_ / core) < second.numerator_ * (first.denominator_ / core);
}

bool operator>=(const Rational& first, const Rational& second) {
  return !(first < second);
}

bool operator<=(const Rational& first, const Rational& second) {
  return (first < second) || (first == second);
}

bool operator>(const Rational& first, const Rational& second) {
  return !(first < second) && (first != second);
}

Rational operator-(const Rational& first, const Rational& second) {
  Rational created = first;
  created -= second;
  return created;
}

Rational operator+(const Rational& first, const Rational& second) {
  Rational created = first - (-second);
  return created;
}

Rational operator*(const Rational& first, const Rational& second) {
  Rational created = first;
  created *= second;
  return created;
}

Rational operator/(const Rational& first, const Rational& second) {
  Rational created = first;
  created /= second;
  return created;
}

std::ostream& operator<<(std::ostream& out, const Rational& rational) {
  return out << rational.toString();
}

std::istream& operator>>(std::istream& in, Rational& rational) {
  BigInteger input;
  in >> input;
  rational = Rational(input);
  return in;
}

// } 'biginteger.h'


struct Element {
  size_t value;
  Element* left = nullptr;
  Element* right = nullptr;

  Element(Element* before, size_t value, Element* after) : value(value), left(before), right(after) {
    before->right = this;
    after->left = this;
  }

  Element(Element* before, size_t value) : value(value), left(before) {
    before->right = this;
  }

  Element(size_t value, Element* after) : value(value), right(after) {
    after->left = this;
  }

  Element(size_t value) : value(value) {}
};

//  Element* operator++(Element*& cur) {
//    cur = cur->right;
//    return cur;
//  }


struct Mask {

  Element* begin = nullptr;
  Element* end = nullptr;

  Mask(ssize_t quantity) {
    begin = new Element(0);
    end = new Element(begin, quantity);
    begin->right = end;

    Element* cur = begin;
    for (int i = 1; i < quantity; i++) {
      cur = new Element(cur, i, cur->right);
    }
  }

  Element* Deactivate(Element* cur) {
    if (cur->left != nullptr) {
      cur->left->right= cur->right;
    } else {
      begin = cur->right;
    }
    cur->right->left = cur->left;
    return cur;
  }

  Element* Activate(Element* cur) {
    if (cur->left != nullptr) {
      cur->left->right = cur;
    } else {
      begin = cur;
    }
    cur->right->left = cur;
    return cur;
  }

  ~Mask() {
    while (begin != end) {
      begin = begin->right;
      delete begin->left;
    }
    delete end;
  }
};


template<size_t N>
class Residue {
 public:
  size_t value;

  Residue(const int key) : value(0) {
    int tool = N;
    value = (key % tool + tool) % tool;
  };

  Residue() : value(0) {};

  explicit Residue(const size_t key) : value(key % N) {};

  explicit operator int() const { return value; };
};

template <size_t N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& it) {
  return out << int(it);
}

template<size_t N, size_t D, bool Stop = (D * D > N)>
struct dividing {
  static const bool result = (N % D == 0) || dividing<N, D + 1>::result;
};

template<size_t N, size_t D>
struct dividing<N, D, true> {
  static const bool result = false;
};

template<size_t N>
struct is_prime {
  static const bool result = !dividing<N, 2>::result;
};

template<>
struct is_prime<0> {
  static const bool result = false;
};

template<size_t N>
Residue<N> operator+(Residue<N> first, Residue<N> second) {
  return Residue<N>((first.value + second.value) % N);
}

template<size_t N>
Residue<N> operator-(Residue<N> first, Residue<N> second) {
  return Residue<N>((N + first.value - second.value) % N);
}

template<size_t N>
Residue<N> operator*(Residue<N> first, Residue<N> second) {
  return Residue<N>((((first.value * second.value) % N) + N) % N);
}

template<size_t N>
size_t gcdex(long long a, long long b, long long &x, long long &y) {  // ax + by = gcd(a, b), b = N
  if (a == 0) {
    x = 0;
    y = 1;
    return x;
  }
  long long x1;
  long long y1;
  gcdex<N>(b % a, a, x1, y1);
  x = y1 - (b / a) * x1;
  y = x1;
  return x;
}

template<size_t N>
Residue<N> Reverse(size_t value) {
  long long x;
  long long y;
  return Residue<N>(gcdex<N>(value, N, x, y) + N);
}

template<size_t N>
Residue<N> operator/(Residue<N> first, Residue<N> second) {
  static_assert(is_prime<N>::result);
  Residue<N> rev = Reverse<N>(second.value);
  return first * rev;
}

template<size_t N>
bool operator==(Residue<N> first, Residue<N> second) {
  return first.value == second.value;
}

template<size_t N>
bool operator!=(Residue<N> first, Residue<N> second) {
  return !(first == second);
}

template<size_t N>
Residue<N>& operator+=(Residue<N>& first, const Residue<N>& second) {
  first = first + second;
  return first;
}

template<size_t N>
Residue<N>& operator-=(Residue<N>& first, const Residue<N>& second) {
  first = first - second;
  return first;
}

template<size_t N>
Residue<N>& operator*=(Residue<N>& first, const Residue<N>& second) {
  first = first * second;
  return first;
}

template<size_t N>
Residue<N>& operator/=(Residue<N>& first, const Residue<N>& second) {
  first = first / second;
  return first;
}


Rational abs(const Rational& it) {
  Rational ret;
  if (it >= 0) {
    ret = it;
  } else {
    ret = -it;
  }
  return ret;
}

template<size_t N>
Residue<N> abs(const Residue<N>&) {
  Residue<N> ret = 0;
  return ret;
}

template<size_t N>
bool operator<(const Residue<N>&, const Residue<N>&) {
  return false;
}

template<size_t N>
bool operator>(const Residue<N>&, const Residue<N>&) {
  return false;
}


template<size_t M, size_t N, typename Field = Rational>
class Matrix {
  std::vector<std::vector<Field>> table;

  Field det(Mask& xs, Mask& ys) const {
    if (xs.begin->right == xs.end) {
      return table[xs.begin->value][ys.begin->value];
    }
    Field su = 0;

    Field sign = 1;
    Element* px = xs.Deactivate(xs.begin);
    Element* py = ys.begin;
    for (; py != ys.end; sign *= Field(-1)) {
      py = ys.Deactivate(py);
      su += table[px->value][py->value] * det(xs, ys) * sign;
      py = ys.Activate(py);
      py = py->right;
    }
    xs.Activate(px);
    return su;
  }

  static void swapRow(std::vector<std::vector<Field>>& mat, size_t first, size_t second) {
    if (first == second) {
      return;
    }
    for (size_t i = 0; i < N; ++i) {
      std::swap(mat[first][i], mat[second][i]);
    }
  }

  static void swapColumn(std::vector<std::vector<Field>>& mat, size_t first, size_t second) {
    if (first == second) {
      return;
    }
    for (size_t i = 0; i < M; ++i) {
      std::swap(mat[i][first], mat[i][second]);
    }
  }

  static void substractRow(std::vector<std::vector<Field>>& mat, size_t object, size_t tool, size_t target) {
    if (mat[tool][target] == Field(0)) {
      return;
    }
    Field coe = mat[object][target] / mat[tool][target];
    for (size_t i = 0; i < N; ++i) {
      mat[object][i] -= mat[tool][i] * coe;
    }
  }

  static void substractColumn(std::vector<std::vector<Field>>& mat, size_t object, size_t tool, size_t target) {
    if (mat[target][tool] == Field(0)) {
      return;
    }
    Field coe = mat[target][object] / mat[target][tool];
    for (size_t i = 0; i < M; ++i) {
      mat[i][object] -= mat[i][tool] * coe;
    }
  }

 public:
  Matrix() : table(std::vector<std::vector<Field>>(M, std::vector<Field>(N, 0))) {}

  Matrix(const Matrix<M, N, Field>& other) : table(std::vector<std::vector<Field>>(M, std::vector<Field>(N, 0))) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        table[i][j] = other[i, j];
      }
    }
  }

  Matrix(std::initializer_list<std::initializer_list<Field>> lst) : table(std::vector<std::vector<Field>>(M, std::vector<Field>(N, 0))) {
    auto row = lst.begin();
    for (size_t i = 0; i < lst.size(); i++, row++) {
      auto cur = row->begin();
      for (size_t j = 0; j < row->size(); j++, cur++) {
        table[i][j] = *cur;
      }
    }
  }

  static Matrix<M, N, Field> unityMatrix() {
    static_assert(M == N);
    Matrix<M, N, Field> created;
    for (size_t i = 0; i < M; ++i) {
      created[i, i] = Field(1);
    }
    return created;
  }

  //  const std::vector<Field>& operator[](size_t index) const {
  //    return const_cast<const std::vector<Field>&>(table[index]);
  //  }
  //
  //  std::vector<Field>& operator[](size_t index) {
  //    return table[index];
  //  }

  void print_int() const {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        std::cerr << table[i][j] << ' ';
      }
      std::cerr << '\n';
    }
    for (size_t i = 0; i < M * N; ++i) {
      std::cerr << Field(i) << ' ';
    }
    std::cerr << '\n';
  }

  Field& operator[](size_t i, size_t j) {
    // print_int();
    return table[i][j];
  }

  const Field& operator[](size_t i, size_t j) const {
    // print_int();
    return const_cast<const Field&>(table[i][j]);;
  }

  Matrix<M, N, Field>& operator=(const Matrix<M, N, Field>& other) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        table[i][j] = other[i, j];
      }
    }
    return *this;
  }

  Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& other) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        *this[i, j] += other[i, j];
      }
    }
    return *this;
  }

  Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& other) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        *this[i, j] -= other[i, j];
      }
    }
    return *this;
  }

  Field det() const {
    // std::cerr << time(0) << " det\n";
    static_assert(M == N);
    Mask xs(M);
    Mask ys(M);
    return det(xs, ys);
    // return LU(xs, ys, N);
  }

  Matrix<N, M, Field> transposed() const {
    auto created = Matrix<N, M, Field>();
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        created[j, i] = table[i][j];
      }
    }
    return created;
  }

  Field LU(Mask& xs, Mask& ys, size_t range) const {
    try {
      std::vector<std::vector<Field>> A;
      for (Element *px = xs.begin; px != xs.end; px = px->right) {
        A.push_back(std::vector<Field>());
        for (Element *py = ys.begin; py != ys.end; py = py->right) {
          A.back().push_back(table[px->value][py->value]);
        }
      }
      std::vector<std::vector<Field>> L = std::vector<std::vector<Field>>(N,
                                                                          std::vector<Field>(N, 0));
      std::vector<std::vector<Field>> U = std::vector<std::vector<Field>>(N,
                                                                          std::vector<Field>(N, 0));
      static_assert(M == N);
      for (size_t i = 0; i < range; i++) {
        for (size_t j = 0; j < range; j++) {
          if (i == j) {
            L[i][j] = 1.0;
          } else {
            L[i][j] = 0.0;
          }
          U[i][j] = 0.0;
        }
      }
      for (size_t k = 0; k < range; k++) {
        U[k][k] = A[k][k];
        for (size_t i = k + 1; i < range; i++) {
          L[i][k] = A[i][k] / U[k][k];
          U[k][i] = A[k][i];
        }
        for (size_t i = k + 1; i < range; i++) {
          for (size_t j = k + 1; j < range; j++) {
            A[i][j] = A[i][j] - L[i][k] * U[k][j];
          }
        }
      }
      Field de = 1;
      for (size_t i = 0; i < range; i++) {
        de *= U[i][i] * L[i][i];
      }
      return de;
    } catch (std::exception& e) {
      // std::cerr << "catch: " << e.what() << '\n';
      return det(xs, ys);
    }
  }

  void swap_rows(std::vector<std::vector<Field>>& A, size_t i, size_t j) const {
    size_t n = A.size();
    for (size_t k = 0; k < n; k++) {
      Field temp = A[i][k];
      A[i][k] = A[j][k];
      A[j][k] = temp;
    }
  }

  std::vector<std::pair<ssize_t, ssize_t>> LU_decomposition(std::vector<std::vector<Field>>& A) const {
    size_t n = A.size();
    std::vector<std::pair<ssize_t, ssize_t>> stack;
    for (size_t j = 0; j < n; j++) {
      size_t pivot = j;
      Field max = abs(A[j][j]);
      for (size_t i = j + 1; i < n; i++) {
        if (abs(A[i][j]) > max) {
          pivot = i;
          max = abs(A[i][j]);
        }
      }
      swap_rows(A, j, pivot);
      if (j != pivot) {
        stack.push_back(std::pair<ssize_t, ssize_t>(j, pivot));
      }
      for (size_t i = j + 1; i < n; i++) {
        A[i][j] /= A[j][j];
        for (size_t k = j + 1; k < n; k++) {
          A[i][k] -= A[i][j] * A[j][k];
        }
      }
    }
    return stack;
  }

  std::vector<Field> forward_substitution(std::vector<std::vector<Field>>& A, std::vector<Field>& b) const {
    size_t n = A.size();
    std::vector<Field> x(n);
    for (size_t i = 0; i < n; i++) {
      Field sum = 0;
      for (size_t j = 0; j < i; j++) {
        sum += A[i][j] * x[j];
      }
      x[i] = b[i] - sum;
    }
    return x;
  }

  std::vector<Field> backward_substitution(std::vector<std::vector<Field>>& A, std::vector<Field>& b) const {
    size_t n = A.size();
    std::vector<Field> x(n);
    for (ssize_t i = n - 1; i >= 0; i--) {
      Field sum = 0;
      for (size_t j = i + 1; j < n; j++) {
        sum += A[i][j] * x[j];
      }
      x[i] = (b[i] - sum) / A[i][i];
    }
    return x;
  }

  std::vector<std::vector<Field>> inverse_matrix(std::vector<std::vector<Field>>& A) const {
    size_t n = A.size();
    std::vector<std::vector<Field>> B(n, std::vector<Field>(n));
    std::vector<std::pair<ssize_t, ssize_t>> stack = LU_decomposition(A);
    for (size_t j = 0; j < n; j++) {
      std::vector<Field> e(n, 0);
      e[j] = 1;
      std::vector<Field> y = forward_substitution(A, e);
      std::vector<Field> x = backward_substitution(A, y);
      for (size_t i = 0; i < n; i++) {
        B[i][j] = x[i];
      }
    }
    for (ssize_t i = stack.size() - 1; i >= 0; --i) {
      for (size_t j = 0; j < N; ++j) {
        std::swap(B[j][stack[i].first], B[j][stack[i].second]);
      }
    }
    return B;
  }

  Matrix<M, N, Field> invertedLU() const {
    static_assert(M == N);
    std::vector<std::vector<Field>> A = table;
    std::vector<std::vector<Field>> B = inverse_matrix(A);
    Matrix<M, N, Field> created;
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        created[i, j] = B[i][j];
      }
    }
    //    Matrix<M, N, Field> correction = *this * created;
    //    Matrix<M, N, Field> real;
    //    for (size_t i = 0; i < M; ++i) {
    //      for (size_t j = 0; j < N; ++j) {
    //        if (correction[i, j] == Field(1)) {
    //          for (size_t k = 0; k < M; ++k) {
    //            real[k, i] = created[k, j];
    //          }
    //          continue;
    //        }
    //      }
    //    }
    return created;
  }

  Matrix<M, N, Field> invertedDet() const {
    static_assert(M == N);
    Mask xs(M);
    Mask ys(M);
    Field coe = det(xs, ys);
    auto created = Matrix<M, N, Field>();

    Field sign = 1;
    Element* px = xs.begin;
    for (; px != xs.end; sign *= Field(-1), sign *= Field(1 - 2 * (M % 2)), px = px->right) {
      xs.Deactivate(px);
      Element* py = ys.begin;
      for (; py != ys.end; sign *= Field(-1), py = py->right) {
        ys.Deactivate(py);
        created[px->value, py->value] = det(xs, ys) * sign / coe;
        ys.Activate(py);
      }
      xs.Activate(px);
    }
    // print_int();
    // created.transposed().print_int();
    return created.transposed();
  }

  Matrix<M, N, Field> inverted() const {
    // std::cerr << time(0) << ' ' << M << 'x' << N << " invert\n";
    static_assert(M == N);
    if (abs(Field(1)) == Field(0)) {
      return invertedDet();
    }
    try {
      return invertedLU();
    } catch (std::exception& e) {
      // std::cerr << e.what() << '\n';
      // std::cerr << time(0) << ' ' << M << 'x' << N << " invert x2\n";
      return invertedDet();
    }
  }

  void invert() {
    static_assert(M == N);
    table = inverted().table;
  }

  const Field trace() const {
    static_assert(M == N);
    Field su = 0;
    for (size_t i = 0; i < M; ++i) {
      su += table[i][i];
    }
    return su;
  }

  size_t rank() const {
    // print_int();
    // std::cerr << time(0) << " rank\n";
    std::vector<std::vector<Field>> buf = table;
    size_t width = std::min(M, N);
    for (size_t r = 0; r < width; ++r) {
      bool found = false;
      size_t x;
      size_t y;
      for (size_t i = r; i < M; ++i) {
        for (size_t j = r; j < N; ++j) {
          if (buf[i][j] != Field(0)) {
            x = i;
            y = j;
            found = true;
            break;
          }
        }
        if (found) {
          break;
        }
      }
      if (!found) {
        return r;
      }
      swapRow(buf, r, x);
      swapColumn(buf, r, y);
      for (size_t i = r + 1; i < M; ++i) {
        substractRow(buf, i, r, r);
      }
      for (size_t i = r + 1; i < N; ++i) {
        substractColumn(buf, i, r, r);
      }
    }
    return width;
  }

  std::array<Field, N> getRow(size_t pos) {
    std::array<Field, N> row;
    for (size_t i = 0; i < N; ++i) {
      row[i] = table[pos][i];
    }
    return row;
  }

  std::array<Field, M> getColumn(size_t pos) {
    std::array<Field, M> column;
    for (size_t i = 0; i < M; ++i) {
      column[i] = table[i][pos];
    }
    return column;
  }
};

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
  auto created = Matrix<M, N, Field>();
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      created[i, j] = first[i, j] + second[i, j];
    }
  }
  return created;
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
  auto created = Matrix<M, N, Field>();
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      created[i, j] = first[i, j] - second[i, j];
    }
  }
  return created;
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator*(const Matrix<M, N, Field>& other, Field coe) {
  auto created = Matrix<M, N, Field>();
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      created[i, j] = other[i, j] * coe;
    }
  }
  return created;
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator*(Field coe, const Matrix<M, N, Field>& other) {
  auto created = Matrix<M, N, Field>();
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      created[i, j] = other[i, j] * coe;
    }
  }
  return created;
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator/(const Matrix<M, N, Field>& other, Field coe) {
  return other * (Field(1) / coe);
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& cur, const Matrix<M, N, Field>& other) {
  cur = cur + other;
  return cur;
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& cur, const Matrix<M, N, Field>& other) {
  cur = cur - other;
  return cur;
}

template<size_t M, size_t N, size_t K, typename Field = Rational>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& first, const Matrix<N, K, Field>& second) {
  // std::cerr << time(0) << " multiply\n";
  auto created = Matrix<M, K, Field>();
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < K; ++j) {
      for (size_t k = 0; k < N; ++k) {
        created[i, j] += first[i, k] * second[k, j];
      }
    }
  }
  return created;
}

template<size_t N, typename Field = Rational>
Matrix<N, N, Field> strassen(const Matrix<N, N, Field>& first, const Matrix<N, N, Field>& second) {
  if (N == 1) {
    return first * second; // base case
  }
  else {
    // split matrices into four submatrices of size N/2
    const size_t half = N / 2;
    Matrix<half, half, Field> A11, A12, A21, A22, B11, B12, B21, B22;
    for (size_t i = 0; i < half; ++i) {
      for (size_t j = 0; j < half; ++j) {
        A11[i, j] = first[i, j];
        A12[i, j] = first[i, j + half];
        A21[i, j] = first[i + half, j];
        A22[i, j] = first[i + half, j + half];
        B11[i, j] = second[i, j];
        B12[i, j] = second[i, j + half];
        B21[i, j] = second[i + half, j];
        B22[i, j] = second[i + half, j + half];
      }
    }
    // compute seven products using recursive calls
    auto P1 = strassen(A11 + A22, B11 + B22);
    auto P2 = strassen(A21 + A22, B11);
    auto P3 = strassen(A11, B12 - B22);
    auto P4 = strassen(A22, B21 - B11);
    auto P5 = strassen(A11 + A12, B22);
    auto P6 = strassen(A21 - A11, B11 + B12);
    auto P7 = strassen(A12 - A22, B21 + B22);
    // combine the results into four submatrices of size N/2
    auto C11 = P1 + P4 - P5 + P7;
    auto C12 = P3 + P5;
    auto C21 = P2 + P4;
    auto C22 = P1 - P2 + P3 + P6;
    // construct the result matrix of size N
    auto result = Matrix<N, N, Field>();
    for (size_t i = 0; i < half; ++i) {
      for (size_t j = 0; j < half; ++j) {
        result[i, j] = C11[i, j];
        result[i, j + half] = C12[i, j];
        result[i + half, j] = C21[i, j];
        result[i + half, j + half] = C22[i, j];
      }
    }
    return result;
  }
}


template<size_t M, size_t N, size_t K, size_t P>
struct check_p {
  static const bool small = P < M || P < N || P < K;
};

template<size_t M, size_t N, size_t K, size_t P, bool check>
struct get_p {
  static const size_t result = 2 * get_p<M, N, K, 2 * P, check_p<M, N, K, 2 * P>::small>::result;
};

template<size_t M, size_t N, size_t K, size_t P>
struct get_p<M, N, K, P, false> {
  static const size_t result = P;
};

template<size_t M, size_t N, typename Field = Rational>
bool operator==(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      if (first[i, j] != second[i, j]) {
        return false;
      }
    }
  }
  return true;
}

template<size_t M, size_t N, typename Field = Rational>
bool operator!=(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
  return !(first == second);
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field>& operator+=(Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
  first = first + second;
  return first;
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field>& operator-=(Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
  first = first - second;
  return first;
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field>& operator*=(Matrix<M, N, Field>& first, const Field& coe) {
  first = first * coe;
  return first;
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator*=(Matrix<M, N, Field>& first, const Matrix<N, N, Field>& second) {
  first = first * second;
  return first;
}


template<size_t M, typename Field = Rational>
class SquareMatrix : public Matrix<M, M, Field> {
 public:
  SquareMatrix() : Matrix<M, M, Field>() {}

  SquareMatrix(std::initializer_list<std::initializer_list<Field>> lst) : Matrix<M, M, Field>(lst) {}
};

template<size_t M, size_t N, typename Field = Rational>
SquareMatrix<M, Field>& operator*=(SquareMatrix<M, Field>& first, const SquareMatrix<M, Field>& second) {
  // std::cerr << "square\n";
  first = first * second;
  return first;
}