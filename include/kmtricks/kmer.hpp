/*****************************************************************************
 *   kmtricks
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

// Packed k-mer representation with alphabet-aware encoding.
// see also kmer_hash.hpp

#pragma once

// std
#include <algorithm>
#include <array>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <vector>
#include <limits>

#include <kmtricks/alphabet.hpp>
#include <kmtricks/exceptions.hpp>

#define DEFAULT_MINIMIZER_KM 1000000000

namespace km
{
constexpr uint8_t MAX_SYMBOL_BITS = 5;
const char bToN[] = {'A', 'C', 'T', 'G'};
const char revN[] = {'T', 'G', 'A', 'C'};
const uint8_t revB[] = {2, 3, 0, 1};
const uint8_t NToB[256] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

const uint8_t rev_table[256] = {
 0xaa, 0xea, 0x2a, 0x6a, 0xba, 0xfa, 0x3a, 0x7a, 0x8a, 0xca, 0xa, 0x4a, 0x9a, 0xda, 0x1a, 0x5a,
 0xae, 0xee, 0x2e, 0x6e, 0xbe, 0xfe, 0x3e, 0x7e, 0x8e, 0xce, 0xe, 0x4e, 0x9e, 0xde, 0x1e, 0x5e,
 0xa2, 0xe2, 0x22, 0x62, 0xb2, 0xf2, 0x32, 0x72, 0x82, 0xc2, 0x2, 0x42, 0x92, 0xd2, 0x12, 0x52,
 0xa6, 0xe6, 0x26, 0x66, 0xb6, 0xf6, 0x36, 0x76, 0x86, 0xc6, 0x6, 0x46, 0x96, 0xd6, 0x16, 0x56,
 0xab, 0xeb, 0x2b, 0x6b, 0xbb, 0xfb, 0x3b, 0x7b, 0x8b, 0xcb, 0xb, 0x4b, 0x9b, 0xdb, 0x1b, 0x5b,
 0xaf, 0xef, 0x2f, 0x6f, 0xbf, 0xff, 0x3f, 0x7f, 0x8f, 0xcf, 0xf, 0x4f, 0x9f, 0xdf, 0x1f, 0x5f,
 0xa3, 0xe3, 0x23, 0x63, 0xb3, 0xf3, 0x33, 0x73, 0x83, 0xc3, 0x3, 0x43, 0x93, 0xd3, 0x13, 0x53,
 0xa7, 0xe7, 0x27, 0x67, 0xb7, 0xf7, 0x37, 0x77, 0x87, 0xc7, 0x7, 0x47, 0x97, 0xd7, 0x17, 0x57,
 0xa8, 0xe8, 0x28, 0x68, 0xb8, 0xf8, 0x38, 0x78, 0x88, 0xc8, 0x8, 0x48, 0x98, 0xd8, 0x18, 0x58,
 0xac, 0xec, 0x2c, 0x6c, 0xbc, 0xfc, 0x3c, 0x7c, 0x8c, 0xcc, 0xc, 0x4c, 0x9c, 0xdc, 0x1c, 0x5c,
 0xa0, 0xe0, 0x20, 0x60, 0xb0, 0xf0, 0x30, 0x70, 0x80, 0xc0, 0x0, 0x40, 0x90, 0xd0, 0x10, 0x50,
 0xa4, 0xe4, 0x24, 0x64, 0xb4, 0xf4, 0x34, 0x74, 0x84, 0xc4, 0x4, 0x44, 0x94, 0xd4, 0x14, 0x54,
 0xa9, 0xe9, 0x29, 0x69, 0xb9, 0xf9, 0x39, 0x79, 0x89, 0xc9, 0x9, 0x49, 0x99, 0xd9, 0x19, 0x59,
 0xad, 0xed, 0x2d, 0x6d, 0xbd, 0xfd, 0x3d, 0x7d, 0x8d, 0xcd, 0xd, 0x4d, 0x9d, 0xdd, 0x1d, 0x5d,
 0xa1, 0xe1, 0x21, 0x61, 0xb1, 0xf1, 0x31, 0x71, 0x81, 0xc1, 0x1, 0x41, 0x91, 0xd1, 0x11, 0x51,
 0xa5, 0xe5, 0x25, 0x65, 0xb5, 0xf5, 0x35, 0x75, 0x85, 0xc5, 0x5, 0x45, 0x95, 0xd5, 0x15, 0x55};

const char bToP[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'};

inline uint8_t protein_to_code(unsigned char c)
{
  switch (c)
  {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'D': case 'd': return 2;
    case 'E': case 'e': return 3;
    case 'F': case 'f': return 4;
    case 'G': case 'g': return 5;
    case 'H': case 'h': return 6;
    case 'I': case 'i': return 7;
    case 'K': case 'k': return 8;
    case 'L': case 'l': return 9;
    case 'M': case 'm': return 10;
    case 'N': case 'n': return 11;
    case 'P': case 'p': return 12;
    case 'Q': case 'q': return 13;
    case 'R': case 'r': return 14;
    case 'S': case 's': return 15;
    case 'T': case 't': return 16;
    case 'V': case 'v': return 17;
    case 'W': case 'w': return 18;
    case 'Y': case 'y': return 19;
    case 'B': case 'b':
    case 'Z': case 'z':
    case 'J': case 'j':
    case 'X': case 'x':
    case '*':
      return 20;
    default:
      return 20;
  }
}

inline uint8_t char_to_code(Alphabet alphabet, unsigned char c)
{
  return alphabet == Alphabet::DNA ? NToB[c] : protein_to_code(c);
}

inline char code_to_char(Alphabet alphabet, uint8_t code)
{
  if (alphabet == Alphabet::DNA)
    return bToN[code & 3];
  if (code >= (sizeof(bToP) / sizeof(bToP[0])))
    return 'X';
  return bToP[code];
}

inline std::string str_rev_comp(const std::string& s, Alphabet alphabet = Alphabet::DNA)
{
  if (alphabet != Alphabet::DNA)
    return s;
  std::string rev;
  for (auto it=s.rbegin(); it!=s.rend(); it++)
    rev.push_back(revN[NToB[*it]]);
  return rev;
}

inline bool is_valid_minimizer(uint64_t value, uint8_t size, Alphabet alphabet)
{
  if (alphabet != Alphabet::DNA)
    return true;
  uint32_t mask1 = std::numeric_limits<uint32_t>::max() >> (((sizeof(uint32_t)*8)-(size*2))+4);
  uint32_t mask01 = static_cast<uint32_t>(0x5555555555555555);
  uint32_t mask00 = mask01 & mask1;
  value = ~((value) | (value >> 2));
  value = ((value >> 1) & value) & mask00;
  return !(value != 0);
}

class Mmer
{
public:
  Mmer() {}
  Mmer(uint64_t value, uint8_t size, Alphabet alphabet = Alphabet::DNA)
  {
    set(value, size, alphabet);
  }

  void set(uint64_t value, uint8_t size, Alphabet alphabet)
  {
    m_size = size;
    m_data = value;
    m_alphabet = alphabet;
  }

  Mmer rev_comp()
  {
    if (m_alphabet != Alphabet::DNA)
      return *this;
    uint64_t rev = 0;
    uint64_t tmp = m_data;
    for (int i=m_size-1; i>=0; i--)
    {
      rev <<= 2;
      rev |= revB[tmp & 3];
      tmp >>= 2;
    }
    return Mmer(rev, m_size, m_alphabet);
  }

  std::string to_string() const
  {
    int i;
    uint64_t tmp = m_data;

    std::string seq(m_size, '\0');
    for (i=m_size-1; i>=0; i--)
    {
        seq[static_cast<size_t>(i)] =
          code_to_char(m_alphabet, tmp & ((1ULL << alphabet_bits(m_alphabet)) - 1));
        tmp = tmp >> alphabet_bits(m_alphabet);
    }
    return seq;
  }

  bool operator>(const Mmer& m) const
  {
    return m_data > m.m_data;
  }

  bool operator<(const Mmer& m) const
  {
    return m_data < m.m_data;
  }

  bool operator==(const Mmer& m) const
  {
    return m_data == m.m_data;
  }

  uint64_t value() const { return m_data; }

private:
  uint64_t m_data {0};
  uint8_t m_size {0};
  Alphabet m_alphabet {Alphabet::DNA};
};

/**
 * @brief Kmer class
 *
 *  This a 2bit representation of a k-mer, using A=0, C=1, T=2, G=3 encoding.
 * It uses uint64_t[] as backend, specializations are provided for MAX_K=32 (uint64_t)
 * and MAX_K=64 (__uint128_t).
 *
 *  @warning Multiple k-mer sizes should not be use at the same time in different k-mer object
 *  with the same specilization because kmer size is maintained as a static inline variable.
 *
 * @tparam MAX_K Maximum k-mer size
 */
template <size_t MAX_K>
class Kmer
{
  using data_ptr64 = const uint64_t*;
  using data_ptr8 = const uint8_t*;

public:
  const static uint16_t m_max_data{(MAX_K * MAX_SYMBOL_BITS + 63) / 64};
  inline static size_t m_kmer_size;
  inline static uint16_t m_n_data;
  inline static Alphabet m_alphabet {Alphabet::DNA};
  inline static uint8_t m_bits_per_symbol {km::alphabet_bits(Alphabet::DNA)};

protected:

  union
  {
    uint64_t m_data[m_max_data];
    uint8_t m_data8[m_max_data * 8];
  };

 public:

  static std::string name()
  {
    std::stringstream ss;
    ss << "Kmer<" << std::to_string(MAX_K) << "> - uint64_t[" << std::to_string(m_max_data) << "]";
    return ss.str();
  }

  static const size_t get_size_bits() { return 8 * sizeof(uint64_t) * m_max_data ;}
  static Alphabet alphabet() { return m_alphabet; }
  static uint8_t bits_per_symbol() { return m_bits_per_symbol; }

  static void set_alphabet(Alphabet alphabet)
  {
    m_alphabet = alphabet;
    m_bits_per_symbol = km::alphabet_bits(alphabet);
    if (m_kmer_size != 0)
    {
      const size_t bits = m_kmer_size * m_bits_per_symbol;
      m_n_data = static_cast<uint16_t>((bits + 63) / 64);
      if (bits > (m_max_data * 64))
        throw ConfigError("k-mer size exceeds storage capacity for current alphabet.");
    }
  }

  static uint32_t slots_for_kmer(size_t kmer_size)
  {
    return static_cast<uint32_t>((kmer_size * m_bits_per_symbol + 63) / 64);
  }

  /***********************
  *    Constructors      *
  ***********************/

  Kmer() { zero(); }
  Kmer(size_t kmer_size) { set_k(kmer_size); }
  Kmer(const std::string& str_kmer) { set_polynom(str_kmer); }

  /***********************
  *    Set data          *
  ***********************/

  void zero()
  {
    std::fill(std::begin(m_data), std::end(m_data), 0);
  }

  void set_k(size_t kmer_size)
  {
    std::fill(std::begin(m_data), std::end(m_data), 0);
    m_kmer_size = kmer_size;
    const size_t bits = m_kmer_size * m_bits_per_symbol;
    m_n_data = static_cast<uint16_t>((bits + 63) / 64);
    if (bits > (m_max_data * 64))
      throw ConfigError("k-mer size exceeds storage capacity for current alphabet.");
  }

  void set64(uint64_t value) { m_data[0] = value; }

  void set64_p(const uint64_t* data)
  {
    for (size_t i=0; i<m_n_data; i++)
      m_data[i] = data[i];
  }

  void set_polynom(const char* data, size_t kmer_size)
  {
    set_k(kmer_size);
    for (size_t i=0; i<kmer_size; i++)
      (*this) = (*this << m_bits_per_symbol) + static_cast<uint64_t>(char_to_code(m_alphabet, data[i]));
  }

  void set_polynom(const std::string& s)
  {
    set_k(s.size());
    for (size_t i=0; i<s.size(); i++)
      (*this) = (*this << m_bits_per_symbol) + static_cast<uint64_t>(char_to_code(m_alphabet, s[i]));
  }

  /***********************
  *    Get data          *
  ***********************/

  uint64_t get64() const { return m_data[0]; }
  data_ptr64 get_data64() const { return m_data; }
  data_ptr8 get_data8() const { return m_data8; }
  uint64_t* get_data64_unsafe() { return m_data; }
#ifdef __SIZEOF_INT128__
  void set128(__uint128_t value)
  {
    m_data[0] = static_cast<uint64_t>(value);
    m_data[1] = static_cast<uint64_t>(value >> 64);
  }
  __uint128_t get128() const
  {
    return (static_cast<__uint128_t>(m_data[1]) << 64) | m_data[0];
  }
#endif

  /***********************
  *    Access data       *
  ***********************/

  uint8_t operator[] (size_t i) const
  {
    const size_t bit_index = i * m_bits_per_symbol;
    const size_t word = bit_index / 64;
    const size_t offset = bit_index % 64;
    uint64_t val = m_data[word] >> offset;
    if (offset + m_bits_per_symbol > 64 && word + 1 < m_max_data)
      val |= m_data[word + 1] << (64 - offset);
    return static_cast<uint8_t>(val & ((1ULL << m_bits_per_symbol) - 1));
  }
  char at(size_t i) const { return code_to_char(m_alphabet, (*this)[m_kmer_size-i-1]); }
  uint8_t at2bit(size_t i) const { return (*this)[m_kmer_size-i-1]; }
  uint8_t byte_at(size_t i) const { return (*this)[m_kmer_size-i-1]; }

  /***********************
  * Comparison operators *
  ***********************/

  bool operator<(const Kmer<MAX_K>& k) const
  {
    for (int i=m_n_data-1; i >= 0; i--)
      if (m_data[i] != k.m_data[i])
        return m_data[i] < k.m_data[i];
    return false;
  }

  bool operator<=(const Kmer<MAX_K>& k) const
  {
    return operator==(k) || operator<(k);
  }

  bool operator>(const Kmer<MAX_K>& k) const
  {
    if (operator==(k)) return false;
    if (operator<(k)) return false;
    return true;
  }

  bool operator>=(const Kmer<MAX_K>& k) const
  {
    return operator==(k) || operator>(k);
  }

  bool operator==(const Kmer<MAX_K>& k) const
  {
    for (uint32_t i = 0; i < m_n_data; i++)
      if (m_data[i] != k.m_data[i])
        return false;
    return true;
  }

  bool operator!=(const Kmer<MAX_K>& k) const
  {
    for (uint32_t i = 0; i < m_n_data; i++)
      if (m_data[i] != k.m_data[i])
        return true;
    return false;
  }

  /***********************
  * Arithmetic operators *
  ***********************/

  Kmer<MAX_K> operator+(const Kmer<MAX_K>& k) const
  {
    Kmer<MAX_K> res;
    uint32_t c = 0;
    for (uint32_t i = 0; i < m_n_data; i++)
    {
      res.m_data[i] = m_data[i] + k.m_data[i] + c;
      c = (res.m_data[i] < m_data[i]) ? 1 : 0;
    }
    return res;
  }

  Kmer<MAX_K> operator+(uint64_t o) const
  {
    Kmer<MAX_K> res;
    uint32_t c = 0;
    res.m_data[0] = m_data[0] + o;
    for (uint32_t i=1; i<m_n_data; i++)
    {
      res.m_data[i] = m_data[i] + c;
      c = (res.m_data[i] < m_data[i]) ? 1 : 0;
    }
    return res;
  }

  Kmer<MAX_K> operator-(const Kmer<MAX_K>& k) const
  {
    Kmer<MAX_K> res;
    uint32_t c = 0;
    for (uint32_t i=0; i<m_max_data; i++)
    {
      res.m_data[i] = m_data[i] - k.m_data[i] - c;
      c = (res.m_data[i] > m_data[i]) ? 1 : 0;
    }
    return res;
  }

  Kmer<MAX_K> operator-(uint64_t o) const
  {
    Kmer<MAX_K> res;
    uint32_t c = 0;
    res.m_data[0] = m_data[0] - o;
    c = (res.m_data[0] > m_data[0]) ? 1 : 0;
    for (uint32_t i=0; i<m_max_data; i++)
    {
      res.m_data[i] = m_data[i] - c;
      c = (res.m_data[i] > m_data[i]) ? 1 : 0;
    }
    return res;
  }

  Kmer<MAX_K> operator*(uint32_t coeff) const
  {
    Kmer<MAX_K> res;
    if (coeff == 0)
      return res;
    if ((coeff & (coeff - 1)) == 0)
    {
      uint32_t shift = 0;
      while ((coeff >> shift) > 1)
        shift++;
      return (*this) << shift;
    }
#ifdef __SIZEOF_INT128__
    uint64_t carry = 0;
    for (size_t i=0; i<m_max_data; i++)
    {
      __uint128_t prod = static_cast<__uint128_t>(m_data[i]) * coeff + carry;
      res.m_data[i] = static_cast<uint64_t>(prod);
      carry = static_cast<uint64_t>(prod >> 64);
    }
#else
    throw ConfigError("Generic k-mer multiplication requires __int128 support.");
#endif
    return res;
  }

  Kmer<MAX_K> operator/(uint32_t coeff) const
  {
    Kmer<MAX_K> res;
    uint64_t r = 0;
    uint32_t mask32 = ~0;
    for (int i=m_max_data-1; i>=0; i--)
    {
      for (int j=1; j>=0; j--)
      {
        uint64_t n = (r << 32) | ((m_data[i] >> (32 * j)) & mask32);
        res.m_data[i] = res.m_data[i] | (((n / coeff) & mask32) << (32 * j));
        r = n % coeff;
      }
    }
    return res;
  }

  uint32_t operator%(uint32_t& coeff) const
  {
    uint64_t r = 0;
    uint32_t mask32 = 0;
    for (int i=m_max_data-1; i>=0; i--)
    {
      for (int j=1; j>=0; j--)
      {
        uint64_t n = (r << 32) | ((m_data[i] >> (32 * j)) & mask32);
        r = n % coeff;
      }
    }
    return static_cast<uint32_t>(r);
  }

  /***********************
  *  Bitwise operators   *
  ***********************/

  Kmer<MAX_K> operator^(const Kmer<MAX_K>& k) const
  {
    Kmer<MAX_K> res;
    for (uint32_t i=0; i<m_max_data; i++)
      res.m_data[i] = m_data[i] ^ k.m_data[i];
    return res;
  }

  Kmer<MAX_K> operator|(const Kmer<MAX_K>& k) const
  {
    Kmer<MAX_K> res;
    for (uint32_t i=0; i<m_max_data; i++)
      res.m_data[i] = m_data[i] | k.m_data[i];
    return res;
  }

  Kmer<MAX_K> operator&(const Kmer<MAX_K>& k) const
  {
    Kmer<MAX_K> res;
    for (uint32_t i=0; i<m_max_data; i++)
      res.m_data[i] = m_data[i] & k.m_data[i];
    return res;
  }

  Kmer<MAX_K> operator&(char o) const
  {
    Kmer<MAX_K> res; res.set64(0);
    res.m_data[0] = m_data[0] & o;
    return res;
  }

  Kmer<MAX_K> operator~() const
  {
    Kmer<MAX_K> res;
    for (uint32_t i=0; i<m_max_data; i++)
      res.m_data[i] = ~m_data[i];
    return res;
  }

  Kmer<MAX_K> operator>>(uint32_t shift) const
  {
    Kmer<MAX_K> res;
    if (shift == 0)
      return *this;
    const size_t wshift = shift / 64;
    const size_t bshift = shift % 64;
    for (size_t i=0; i<m_max_data; i++)
    {
      if (i + wshift >= m_max_data)
      {
        res.m_data[i] = 0;
        continue;
      }
      uint64_t val = m_data[i + wshift] >> bshift;
      if (bshift && i + wshift + 1 < m_max_data)
        val |= m_data[i + wshift + 1] << (64 - bshift);
      res.m_data[i] = val;
    }
    return res;
  }

  Kmer<MAX_K> operator<<(uint32_t shift) const
  {
    Kmer<MAX_K> res;
    if (shift == 0)
      return *this;
    const size_t wshift = shift / 64;
    const size_t bshift = shift % 64;
    for (size_t i=m_max_data; i-- > 0; )
    {
      if (i < wshift)
      {
        res.m_data[i] = 0;
        continue;
      }
      uint64_t val = m_data[i - wshift] << bshift;
      if (bshift && i > wshift)
        val |= m_data[i - wshift - 1] >> (64 - bshift);
      res.m_data[i] = val;
    }
    return res;
  }

  /***********************
  * Assignment operators *
  ***********************/

  Kmer<MAX_K>& operator+=(const Kmer<MAX_K>& k) { *this = *this + k; return *this; }
  Kmer<MAX_K>& operator-=(const Kmer<MAX_K>& k) { *this = *this - k; return *this; }

  Kmer<MAX_K>& operator*=(uint32_t coeff) { *this = *this * coeff; return *this; }
  Kmer<MAX_K>& operator/=(uint32_t coeff) { *this = *this / coeff; return *this; }

  Kmer<MAX_K>& operator&=(const Kmer<MAX_K>& k)
  {
    for (uint32_t i=0; i<m_max_data; i++)
      m_data[i] &= k.m_data[i];
    return *this;
  }

  Kmer<MAX_K>& operator|=(const Kmer<MAX_K>& k)
  {
    for (uint32_t i=0; i<m_max_data; i++)
      m_data[i] |= k.m_data[i];
    return *this;
  }

  Kmer<MAX_K>& operator^=(const Kmer<MAX_K>& k)
  {
    for (uint32_t i=0; i<m_max_data; i++)
      m_data[i] ^= k.m_data[i];
    return *this;
  }

  Kmer<MAX_K>& operator<<=(uint32_t shift) { *this = *this << shift; return *this; }
  Kmer<MAX_K>& operator>>=(uint32_t shift) { *this = *this >> shift; return *this; }

  /***********************
  *   Kmer operations    *
  ***********************/

  Kmer<MAX_K> rev_comp() const
  {
    if (m_alphabet != Alphabet::DNA)
      return *this;
    Kmer<MAX_K> kmer; kmer.set_k(m_kmer_size);
    for (size_t i=0; i<8*m_n_data; i++)
    {
      kmer.m_data8[8*m_n_data-1-i] = rev_table[m_data8[i]];
    }
    const size_t symbols_per_word = 64 / m_bits_per_symbol;
    return (kmer >> (m_bits_per_symbol * (symbols_per_word * m_n_data - m_kmer_size)));
  }

  Kmer<MAX_K> canonical() const
  {
    if (m_alphabet != Alphabet::DNA)
      return *this;
    Kmer<MAX_K> kmer = rev_comp();
    return (kmer < *this) ? kmer : *this;
  }

  /***************************
  *   Text representation    *
  ***************************/

  std::string to_string() const
  {
    char seq[m_kmer_size + 1];
    for (size_t i=0; i<m_kmer_size; i++)
    {
      seq[m_kmer_size-i-1] = code_to_char(m_alphabet, (*this)[i]);
    }
    seq[m_kmer_size] = '\0';
    return seq;
  }

  std::string to_bit_string() const
  {
    std::stringstream ss;
    for (uint16_t i = 0; i < m_n_data; i++)
      ss << i << " " << std::bitset<sizeof(uint64_t) * 8>(m_data[i]).to_string() << "\n";
    return ss.str();
  }

  /***************************
  *   Stream methods         *
  ***************************/

  void dump(std::ostream& stream)
  {
    stream.write(reinterpret_cast<char*>(m_data), m_n_data * sizeof(uint64_t));
  }

  void load(std::istream& stream)
  {
    stream.read(reinterpret_cast<char*>(m_data), m_n_data * sizeof(uint64_t));
  }

  std::vector<Mmer> mmers(uint8_t size) const
  {
    const size_t nb_mmers = m_kmer_size - size + 1;
    std::vector<Mmer> mmers(nb_mmers);
    for (size_t i=0; i<nb_mmers; i++)
    {
      uint64_t value = 0;
      for (size_t j=i; j<i+size; j++)
      {
        value <<= m_bits_per_symbol;
        value |= byte_at(j);
      }
      mmers[i].set(value, size, m_alphabet);
    }
    return mmers;
  }

  Mmer minimizer(uint8_t size)
  {
    const uint32_t max_bits = size * m_bits_per_symbol;
    if (max_bits >= 64)
      throw ConfigError("Minimizer size exceeds supported encoding width.");
    uint64_t def = ((uint64_t)1 << max_bits) - 1;
    const size_t nb_mmers = m_kmer_size - size + 1;
    Mmer minim(std::numeric_limits<uint64_t>::max(), size, m_alphabet);
    for (size_t i=0; i<nb_mmers; i++)
    {
      uint64_t value = 0;
      for (size_t j=i; j<i+size; j++)
      {
        value <<= m_bits_per_symbol;
        value |= byte_at(j);
      }
      uint64_t rev = value;
      if (m_alphabet == Alphabet::DNA)
      {
        rev = 0;
        uint64_t tmp = value;
        for (int j=size-1; j>=0; j--)
        {
          rev <<= 2;
          rev |= revB[tmp & 3];
          tmp >>= 2;
        }
      }
      Mmer tmp_minim(rev < value ? rev : value, size, m_alphabet);
      if (is_valid_minimizer(tmp_minim.value(), size, m_alphabet))
      {
        if (tmp_minim < minim)
          minim = tmp_minim;
      }
      else
      {
        tmp_minim = Mmer(def, size, m_alphabet);
        if (tmp_minim < minim)
        {
          minim = tmp_minim;
        }
      }
    }
    return minim;
  }
};

/**
 * @brief Kmer class with a count value
 *
 * @tparam MAX_K
 * @tparam count_type
 */
template<size_t MAX_K, typename count_type>
class CKmer : public Kmer<MAX_K>
{
public:
  CKmer() {}

  CKmer(const std::string& s, count_type count)
    : Kmer<MAX_K>(s), m_count(count)
  {}

  void set_count(count_type count)
  {
    m_count = count;
  }

  count_type get_count() const
  {
    return m_count;
  }

  void dump_with_count(std::ostream& stream)
  {
    this->dump(stream);
    stream.write(reinterpret_cast<char*>(&m_count), sizeof(m_count));
  }

  void load_with_count(std::istream& stream)
  {
    this->load(stream);
    stream.read(reinterpret_cast<char*>(&m_count), sizeof(m_count));
  }

private:
  count_type m_count;
};

};  // namespace km
