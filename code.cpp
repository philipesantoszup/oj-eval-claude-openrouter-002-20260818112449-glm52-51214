#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

// Do not use "using namespace std;"

namespace sjtu {
class int2048 {
  // Storage: little-endian limbs in base 1e9, plus a sign flag.
  // sign == false  -> value >= 0
  // sign == true   -> value < 0
  // The value 0 is always stored with sign == false and an empty limb vector.
private:
  static const unsigned int BASE = 1000000000u; // 1e9
  static const int WIDTH = 9;
  bool sign;
  std::vector<unsigned int> d; // limbs, little-endian; empty == 0

  // ---- internal helpers ----
  void trim();                       // remove leading zero limbs, fix sign for 0
  bool is_zero() const;              // true iff value == 0
  int cmp_abs(const int2048 &) const; // -1/0/1 comparing |*this| vs |rhs|
  static int2048 add_abs(const int2048 &, const int2048 &); // |a|+|b|, sign+
  static int2048 sub_abs(const int2048 &, const int2048 &); // |a|-|b|, a>=b
  static int2048 mul_abs(const int2048 &, const int2048 &); // |a|*|b|, sign+
  // magnitude division: q = |a|/|b|, r = |a|%|b|, 0 <= r < |b|
  static void divmod_abs(const int2048 &a, const int2048 &b,
                         int2048 &q, int2048 &r);

public:
  // Constructors
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);
  // move constructor / assignment (added for efficiency, compatible)
  int2048(int2048 &&) noexcept;
  int2048 &operator=(int2048 &&) noexcept;

  // The parameter types of the following functions are for reference only, you can choose to use constant references or not
  // If needed, you can add other required functions yourself
  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string &);
  // Output the stored big integer, no need for newline
  void print();

  // Add a big integer
  int2048 &add(const int2048 &);
  // Return the sum of two big integers
  friend int2048 add(int2048, const int2048 &);

  // Subtract a big integer
  int2048 &minus(const int2048 &);
  // Return the difference of two big integers
  friend int2048 minus(int2048, const int2048 &);

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);

  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);

  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);

  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);

  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);

  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);

  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);

  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};
} // namespace sjtu

#endif


namespace sjtu {

// ============================================================
// internal NTT (3-prime CRT) for large multiplication
// ============================================================
namespace {

typedef unsigned long long u64;
typedef unsigned int u32;
typedef __int128 i128;
typedef unsigned __int128 u128;

const u64 MOD0 = 998244353u;  // 119*2^23+1, g=3
const u64 MOD1 = 1004535809u; // 479*2^21+1, g=3
const u64 MOD2 = 469762049u;  // 7*2^26+1, g=3

u64 mpow(u64 a, u64 e, u64 mod) {
  u64 r = 1 % mod;
  a %= mod;
  while (e) {
    if (e & 1u) r = r * a % mod;
    a = a * a % mod;
    e >>= 1;
  }
  return r;
}

// In-place NTT over a prime with primitive root g=3.
// a.size() must be a power of two. invert=true => inverse transform.
void ntt(u32 *a, int n, u64 mod, bool invert) {
  // bit-reversal permutation
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      u32 t = a[i];
      a[i] = a[j];
      a[j] = t;
    }
  }
  for (int len = 2; len <= n; len <<= 1) {
    u64 w = invert ? mpow(3, mod - 1 - (mod - 1) / len, mod)
                   : mpow(3, (mod - 1) / len, mod);
    int half = len >> 1;
    for (int i = 0; i < n; i += len) {
      u64 wn = 1;
      u32 *p = a + i;
      u32 *q = a + i + half;
      for (int j = 0; j < half; ++j) {
        u64 u = p[j];
        u64 v = (u64)q[j] * wn % mod;
        u64 s = u + v;
        if (s >= mod) s -= mod;
        u64 t = u + mod - v;
        if (t >= mod) t -= mod;
        p[j] = (u32)s;
        q[j] = (u32)t;
        wn = wn * w % mod;
      }
    }
  }
  if (invert) {
    u64 inv_n = mpow((u64)n, mod - 2, mod);
    for (int i = 0; i < n; ++i) a[i] = (u32)((u64)a[i] * inv_n % mod);
  }
}

// CRT constants (precomputed at compile time)
// M01 = MOD0*MOD1, inv0_1 = MOD0^(MOD1-2) mod MOD1, inv01_2 = M01^(MOD2-2) mod MOD2
const u64 M01 = MOD0 * MOD1;
const u64 INV0_MOD1 = mpow(MOD0 % MOD1, MOD1 - 2, MOD1);
const u64 INV01_MOD2 = mpow(M01 % MOD2, MOD2 - 2, MOD2);
const u64 M01_MOD2 = M01 % MOD2;
const u64 MOD0_MOD2 = MOD0 % MOD2;

// reconstruct full value from three residues (< MOD0*MOD1*MOD2 ~ 2^88)
inline u128 crt(u32 r0, u32 r1, u32 r2) {
  // v01 = r0 + MOD0 * t1, where t1 = (r1 - r0) * INV0_MOD1 mod MOD1
  u64 t1 = ((u64)r1 + MOD1 - (u64)(r0 % MOD1)) % MOD1;
  t1 = t1 * INV0_MOD1 % MOD1;
  u64 v01 = (u64)r0 + MOD0 * t1; // < MOD0*MOD1 = M01
  // v = v01 + M01 * t2, t2 = (r2 - v01) * INV01_MOD2 mod MOD2
  u64 v01_mod2 = v01 % MOD2;
  u64 t2 = ((u64)r2 + MOD2 - v01_mod2) % MOD2;
  t2 = t2 * INV01_MOD2 % MOD2;
  return (u128)v01 + (u128)M01 * (u128)t2;
}

} // namespace

// ============================================================
// basic helpers
// ============================================================

void int2048::trim() {
  while (!d.empty() && d.back() == 0u) d.pop_back();
  if (d.empty()) sign = false;
}

bool int2048::is_zero() const { return d.empty(); }

int int2048::cmp_abs(const int2048 &rhs) const {
  int la = (int)d.size(), lb = (int)rhs.d.size();
  if (la != lb) return la < lb ? -1 : 1;
  for (int i = la - 1; i >= 0; --i) {
    if (d[i] != rhs.d[i]) return d[i] < rhs.d[i] ? -1 : 1;
  }
  return 0;
}

// magnitude addition
int2048 int2048::add_abs(const int2048 &a, const int2048 &b) {
  int2048 r;
  r.sign = false;
  int la = (int)a.d.size(), lb = (int)b.d.size();
  int n = la > lb ? la : lb;
  r.d.reserve(n + 1);
  u64 carry = 0;
  for (int i = 0; i < n; ++i) {
    u64 cur = carry;
    if (i < la) cur += a.d[i];
    if (i < lb) cur += b.d[i];
    r.d.push_back((u32)(cur % BASE));
    carry = cur / BASE;
  }
  if (carry) r.d.push_back((u32)carry);
  r.trim();
  return r;
}

// magnitude subtraction, assumes |a| >= |b|
int2048 int2048::sub_abs(const int2048 &a, const int2048 &b) {
  int2048 r;
  r.sign = false;
  int la = (int)a.d.size(), lb = (int)b.d.size();
  r.d.reserve(la);
  long long borrow = 0;
  for (int i = 0; i < la; ++i) {
    long long cur = (long long)a.d[i] - borrow;
    if (i < lb) cur -= (long long)b.d[i];
    if (cur < 0) {
      cur += BASE;
      borrow = 1;
    } else {
      borrow = 0;
    }
    r.d.push_back((u32)cur);
  }
  r.trim();
  return r;
}

// ============================================================
// constructors / assignment / IO
// ============================================================

int2048::int2048() : sign(false) {}

int2048::int2048(long long x) : sign(false) {
  if (x < 0) {
    sign = true;
    // careful: -(min long long) overflows; use unsigned
    u64 ux = (u64)(-(x + 1)) + 1u; // safe absolute value
    while (ux) {
      d.push_back((u32)(ux % BASE));
      ux /= BASE;
    }
  } else if (x > 0) {
    u64 ux = (u64)x;
    while (ux) {
      d.push_back((u32)(ux % BASE));
      ux /= BASE;
    }
  }
}

int2048::int2048(const std::string &s) : sign(false) { read(s); }

int2048::int2048(const int2048 &o) : sign(o.sign), d(o.d) {}

int2048::int2048(int2048 &&o) noexcept : sign(o.sign), d(std::move(o.d)) {
  o.sign = false;
}

int2048 &int2048::operator=(const int2048 &o) {
  if (this != &o) {
    sign = o.sign;
    d = o.d;
  }
  return *this;
}

int2048 &int2048::operator=(int2048 &&o) noexcept {
  if (this != &o) {
    sign = o.sign;
    d = std::move(o.d);
    o.sign = false;
  }
  return *this;
}

void int2048::read(const std::string &s) {
  d.clear();
  sign = false;
  int pos = 0, len = (int)s.size();
  if (len == 0) return;
  if (s[0] == '-') {
    sign = true;
    pos = 1;
  } else if (s[0] == '+') {
    pos = 1;
  }
  // skip leading zeros
  while (pos < len && s[pos] == '0') ++pos;
  if (pos >= len) {
    // all zeros (or empty after sign)
    sign = false;
    return;
  }
  // parse from the most significant end in chunks of WIDTH
  int end = len;
  int start = pos;
  // number of digits = end - start
  // process full 9-digit chunks from the right
  int firstLen = (end - start) % WIDTH;
  int idx = start;
  if (firstLen > 0) {
    u32 chunk = 0;
    for (int k = 0; k < firstLen; ++k) {
      chunk = chunk * 10u + (u32)(s[idx + k] - '0');
    }
    d.push_back(chunk);
    idx += firstLen;
  }
  while (idx < end) {
    u32 chunk = 0;
    for (int k = 0; k < WIDTH; ++k) {
      chunk = chunk * 10u + (u32)(s[idx + k] - '0');
    }
    d.push_back(chunk);
    idx += WIDTH;
  }
  // d was built from most-significant to least-significant; reverse to little-endian
  int i = 0, j = (int)d.size() - 1;
  while (i < j) {
    u32 t = d[i];
    d[i] = d[j];
    d[j] = t;
    ++i;
    --j;
  }
  trim();
}

void int2048::print() {
  if (d.empty()) {
    std::putchar('0');
    return;
  }
  if (sign) std::putchar('-');
  // most significant limb without leading zeros
  std::printf("%u", d.back());
  for (int i = (int)d.size() - 2; i >= 0; --i) {
    std::printf("%09u", d[i]);
  }
}

// ============================================================
// signed add / minus
// ============================================================

int2048 &int2048::add(const int2048 &rhs) {
  if (sign == rhs.sign) {
    bool s = sign;
    *this = add_abs(*this, rhs);
    sign = s; // add_abs returns sign=false; restore original sign
  } else {
    int c = cmp_abs(rhs);
    if (c == 0) {
      d.clear();
      sign = false;
    } else if (c > 0) {
      bool s = sign;
      *this = sub_abs(*this, rhs);
      sign = s;
    } else {
      bool s = rhs.sign;
      *this = sub_abs(rhs, *this);
      sign = s;
    }
  }
  trim();
  return *this;
}

int2048 &int2048::minus(const int2048 &rhs) {
  // a - b == a + (-b); flip rhs sign conceptually
  if (sign != rhs.sign) {
    // different signs: magnitudes add, result keeps *this sign
    bool s = sign;
    *this = add_abs(*this, rhs);
    sign = s;
  } else {
    // same sign: subtract magnitudes
    int c = cmp_abs(rhs);
    if (c == 0) {
      d.clear();
      sign = false;
    } else if (c > 0) {
      bool s = sign;
      *this = sub_abs(*this, rhs);
      sign = s;
    } else {
      bool s = !rhs.sign; // result takes opposite of rhs sign
      *this = sub_abs(rhs, *this);
      sign = s;
    }
  }
  trim();
  return *this;
}

// ============================================================
// multiplication (schoolbook + 3-prime NTT)
// ============================================================

int2048 int2048::mul_abs(const int2048 &a, const int2048 &b) {
  if (a.is_zero() || b.is_zero()) return int2048();
  int na = (int)a.d.size(), nb = (int)b.d.size();
  int2048 r;
  r.sign = false;
  long long prod = (long long)na * (long long)nb;
  if (prod <= 100000) {
    // schoolbook
    r.d.assign(na + nb, 0u);
    for (int i = 0; i < na; ++i) {
      u64 carry = 0;
      u64 ai = a.d[i];
      int base = i;
      for (int j = 0; j < nb; ++j) {
        u64 cur = (u64)r.d[base + j] + ai * b.d[j] + carry;
        r.d[base + j] = (u32)(cur % BASE);
        carry = cur / BASE;
      }
      int k = base + nb;
      while (carry) {
        u64 cur = (u64)r.d[k] + carry;
        r.d[k] = (u32)(cur % BASE);
        carry = cur / BASE;
        ++k;
      }
    }
    r.trim();
    return r;
  }
  // NTT path (3-prime CRT)
  int res_len = na + nb; // result has at most res_len limbs
  int L = 1;
  while (L < res_len) L <<= 1;
  std::vector<u32> A0(L, 0), A1(L, 0), A2(L, 0);
  std::vector<u32> B0(L, 0), B1(L, 0), B2(L, 0);
  for (int i = 0; i < na; ++i) {
    u32 v = a.d[i];
    A0[i] = (u32)((u64)v % MOD0);
    A1[i] = (u32)((u64)v % MOD1);
    A2[i] = (u32)((u64)v % MOD2);
  }
  for (int i = 0; i < nb; ++i) {
    u32 v = b.d[i];
    B0[i] = (u32)((u64)v % MOD0);
    B1[i] = (u32)((u64)v % MOD1);
    B2[i] = (u32)((u64)v % MOD2);
  }
  ntt(A0.data(), L, MOD0, false);
  ntt(B0.data(), L, MOD0, false);
  ntt(A1.data(), L, MOD1, false);
  ntt(B1.data(), L, MOD1, false);
  ntt(A2.data(), L, MOD2, false);
  ntt(B2.data(), L, MOD2, false);
  for (int i = 0; i < L; ++i) {
    A0[i] = (u32)((u64)A0[i] * B0[i] % MOD0);
    A1[i] = (u32)((u64)A1[i] * B1[i] % MOD1);
    A2[i] = (u32)((u64)A2[i] * B2[i] % MOD2);
  }
  ntt(A0.data(), L, MOD0, true);
  ntt(A1.data(), L, MOD1, true);
  ntt(A2.data(), L, MOD2, true);
  // CRT + carry propagation
  r.d.clear();
  r.d.reserve(res_len);
  u128 carry = 0;
  for (int i = 0; i < res_len; ++i) {
    u128 v = crt(A0[i], A1[i], A2[i]) + carry;
    r.d.push_back((u32)(v % BASE));
    carry = v / BASE;
  }
  while (carry) {
    r.d.push_back((u32)(carry % BASE));
    carry /= BASE;
  }
  r.trim();
  return r;
}

int2048 &int2048::operator*=(const int2048 &rhs) {
  bool s = sign ^ rhs.sign;
  int2048 res = mul_abs(*this, rhs);
  *this = std::move(res);
  sign = s;
  trim();
  return *this;
}

// ============================================================
// division (Knuth Algorithm D) + floor modulo
// ============================================================

void int2048::divmod_abs(const int2048 &A, const int2048 &B,
                         int2048 &Q, int2048 &R) {
  Q.sign = false;
  R.sign = false;
  if (B.is_zero()) { // undefined by spec; guard against crashes
    Q.d.clear();
    R.d.clear();
    return;
  }
  int c = A.cmp_abs(B);
  if (c < 0) {
    Q.d.clear();
    R.d = A.d; // |A|
    R.sign = false;
    R.trim();
    return;
  }
  if (c == 0) {
    Q.d.clear();
    Q.d.push_back(1u);
    R.d.clear();
    return;
  }
  int n = (int)B.d.size();
  int m = (int)A.d.size() - n; // m >= 0
  if (n == 1) {
    u32 dv = B.d[0];
    Q.d.assign(A.d.size(), 0u);
    u64 rem = 0;
    for (int i = (int)A.d.size() - 1; i >= 0; --i) {
      u64 cur = rem * BASE + A.d[i];
      Q.d[i] = (u32)(cur / dv);
      rem = cur % dv;
    }
    R.d.clear();
    if (rem) R.d.push_back((u32)rem);
    Q.trim();
    R.trim();
    return;
  }
  // normalize: d = BASE / (B.d[n-1] + 1)
  u32 d = (u32)(BASE / ((u64)B.d[n - 1] + 1));
  std::vector<u32> un(m + n + 1, 0u), vn(n, 0u);
  // vn = B * d
  {
    u64 carry = 0;
    for (int i = 0; i < n; ++i) {
      u64 cur = (u64)B.d[i] * d + carry;
      vn[i] = (u32)(cur % BASE);
      carry = cur / BASE;
    }
    // carry must be 0 here (guaranteed by choice of d)
  }
  // un = A * d (length m+n+1)
  {
    u64 carry = 0;
    int la = (int)A.d.size();
    for (int i = 0; i < la; ++i) {
      u64 cur = (u64)A.d[i] * d + carry;
      un[i] = (u32)(cur % BASE);
      carry = cur / BASE;
    }
    un[la] = (u32)carry;
  }
  Q.d.assign(m + 1, 0u);
  for (int j = m; j >= 0; --j) {
    u64 num = (u64)un[j + n] * BASE + un[j + n - 1];
    u64 qhat = num / vn[n - 1];
    u64 rhat = num % vn[n - 1];
    while (qhat >= BASE ||
           qhat * vn[n - 2] > rhat * BASE + un[j + n - 2]) {
      --qhat;
      rhat += vn[n - 1];
      if (rhat >= BASE) break;
    }
    // multiply vn by qhat and subtract from un[j..j+n]
    long long borrow = 0;
    u64 carry = 0;
    for (int i = 0; i < n; ++i) {
      u64 p = qhat * vn[i] + carry;
      carry = p / BASE;
      u64 pl = p % BASE;
      long long sub = (long long)un[j + i] - (long long)pl - borrow;
      if (sub < 0) {
        sub += (long long)BASE;
        borrow = 1;
      } else {
        borrow = 0;
      }
      un[j + i] = (u32)sub;
    }
    long long sub = (long long)un[j + n] - (long long)carry - borrow;
    if (sub < 0) {
      sub += (long long)BASE;
      borrow = 1;
    } else {
      borrow = 0;
    }
    un[j + n] = (u32)sub;
    if (borrow) {
      --qhat;
      u64 carry2 = 0;
      for (int i = 0; i < n; ++i) {
        u64 cur = (u64)un[j + i] + vn[i] + carry2;
        un[j + i] = (u32)(cur % BASE);
        carry2 = cur / BASE;
      }
      u64 top = (u64)un[j + n] + carry2;
      un[j + n] = (u32)(top >= BASE ? top - BASE : top);
    }
    Q.d[j] = (u32)qhat;
  }
  // remainder = un[0..n-1] / d
  R.d.assign(n, 0u);
  {
    u64 rem = 0;
    for (int i = n - 1; i >= 0; --i) {
      u64 cur = rem * BASE + un[i];
      R.d[i] = (u32)(cur / d);
      rem = cur % d;
    }
  }
  Q.trim();
  R.trim();
}

int2048 &int2048::operator/=(const int2048 &rhs) {
  int2048 q, r;
  divmod_abs(*this, rhs, q, r);
  bool sa = sign, sb = rhs.sign;
  bool qsign;
  if (sa == sb) {
    qsign = false; // floor quotient is non-negative
  } else {
    if (r.is_zero()) {
      qsign = q.is_zero() ? false : true;
    } else {
      int2048 one((long long)1);
      q = add_abs(q, one); // qmag + 1
      qsign = true;
    }
  }
  *this = std::move(q);
  sign = qsign;
  trim();
  return *this;
}

int2048 &int2048::operator%=(const int2048 &rhs) {
  int2048 q, r;
  divmod_abs(*this, rhs, q, r);
  bool sa = sign, sb = rhs.sign;
  bool rzero = r.is_zero();
  if (sa == sb) {
    if (rzero) {
      d.clear();
      sign = false;
    } else {
      *this = std::move(r);
      sign = sb;
    }
  } else {
    if (rzero) {
      d.clear();
      sign = false;
    } else {
      int2048 babs;
      babs.sign = false;
      babs.d = rhs.d; // |B|
      *this = sub_abs(babs, r); // |B| - rmag
      sign = sb;
    }
  }
  trim();
  return *this;
}

// ============================================================
// unary, compound, comparison, stream operators
// ============================================================

int2048 int2048::operator+() const { return *this; }

int2048 int2048::operator-() const {
  int2048 r(*this);
  if (!r.d.empty()) r.sign = !r.sign;
  return r;
}

int2048 &int2048::operator+=(const int2048 &rhs) { return add(rhs); }
int2048 &int2048::operator-=(const int2048 &rhs) { return minus(rhs); }

int2048 add(int2048 a, const int2048 &b) { a.add(b); return a; }
int2048 minus(int2048 a, const int2048 &b) { a.minus(b); return a; }

int2048 operator+(int2048 a, const int2048 &b) { a += b; return a; }
int2048 operator-(int2048 a, const int2048 &b) { a -= b; return a; }
int2048 operator*(int2048 a, const int2048 &b) { a *= b; return a; }
int2048 operator/(int2048 a, const int2048 &b) { a /= b; return a; }
int2048 operator%(int2048 a, const int2048 &b) { a %= b; return a; }

bool operator==(const int2048 &a, const int2048 &b) {
  if (a.sign != b.sign) return false; // 0 stored as sign=false on both sides
  return a.d == b.d;
}
bool operator!=(const int2048 &a, const int2048 &b) { return !(a == b); }
bool operator<(const int2048 &a, const int2048 &b) {
  if (a.sign != b.sign) return a.sign; // neg < pos
  int c = a.cmp_abs(b);
  if (!a.sign) return c < 0; // both >= 0
  return c > 0;              // both < 0: larger magnitude is smaller
}
bool operator>(const int2048 &a, const int2048 &b) { return b < a; }
bool operator<=(const int2048 &a, const int2048 &b) { return !(a > b); }
bool operator>=(const int2048 &a, const int2048 &b) { return !(a < b); }

namespace {
inline void append_u32(std::string &out, u32 v) {
  if (v == 0) {
    out.push_back('0');
    return;
  }
  char buf[10];
  int n = 0;
  while (v) {
    buf[n++] = (char)('0' + (int)(v % 10u));
    v /= 10u;
  }
  while (n) out.push_back(buf[--n]);
}
inline void append_padded9(std::string &out, u32 v) {
  char buf[9];
  for (int i = 8; i >= 0; --i) {
    buf[i] = (char)('0' + (int)(v % 10u));
    v /= 10u;
  }
  out.append(buf, 9);
}
} // namespace

std::istream &operator>>(std::istream &is, int2048 &x) {
  std::string s;
  is >> s;
  x.read(s);
  return is;
}

std::ostream &operator<<(std::ostream &os, const int2048 &x) {
  if (x.d.empty()) {
    os.put('0');
    return os;
  }
  std::string out;
  out.reserve(x.d.size() * 9 + 2);
  if (x.sign) out.push_back('-');
  append_u32(out, x.d.back());
  for (int i = (int)x.d.size() - 2; i >= 0; --i) {
    append_padded9(out, x.d[i]);
  }
  os << out;
  return os;
}

} // namespace sjtu
