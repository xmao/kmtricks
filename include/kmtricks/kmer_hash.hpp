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

#pragma once
#include <kmtricks/kmer.hpp>

#ifdef WITH_XXHASH
#include <xxhash.h>
#endif

namespace km {

template<size_t MAX_K>
struct IKHasher
{
  virtual uint64_t operator()(const Kmer<MAX_K>& kmer, uint64_t seed = 0) const = 0;
};

template<int HashFunction>
struct KmerHashers {};

template<>
struct KmerHashers<0>
{
  static std::string name() { return "KmerHashers<0> - Folly hash"; }
  template<size_t MAX_K>

  struct Hasher : public IKHasher<MAX_K>
  {
    static std::string name()
    {
      return "KmerHashers<0>::Hasher<MAX_K=" + std::to_string(MAX_K) + ">";
    }

    uint64_t operator()(const Kmer<MAX_K>& kmer, uint64_t seed = 0) const final
    {
      uint64_t _hash = seed, _key = 0;
      for (size_t i=0; i<kmer.m_n_data; i++)
      {
        _key = kmer.get_data64()[i];
        _key = (~_key) + (_key << 21);
        _key = _key ^ (_key >> 24);
        _key = _key + (_key << 3) + (_key << 8);
        _key = _key ^ (_key >> 14);
        _key = _key + (_key << 2) + (_key << 4);
        _key = _key ^ (_key >> 28);
        _key = _key + (_key << 31);
        _hash ^= _key;
      }
      return _hash;
    }
  };

  template<size_t MAX_K>
  struct WinHasher : public IKHasher<MAX_K>
  {
    WinHasher(uint64_t p, uint64_t w) : p(p), w(w) {}
    static std::string name()
    {
      return "KmerHashers<0>::WinHasher<MAX_K=" + std::to_string(MAX_K) + ">";
    }

    uint64_t operator()(const Kmer<MAX_K>& kmer, uint64_t seed = 0) const final
    {
      uint64_t _hash = seed, _key = 0;
      for (size_t i=0; i<kmer.m_n_data; i++)
      {
        _key = kmer.get_data64()[i];
        _key = (~_key) + (_key << 21);
        _key = _key ^ (_key >> 24);
        _key = _key + (_key << 3) + (_key << 8);
        _key = _key ^ (_key >> 14);
        _key = _key + (_key << 2) + (_key << 4);
        _key = _key ^ (_key >> 28);
        _key = _key + (_key << 31);
        _hash ^= _key;
      }
      return (_hash % w) + (w * p);
    }
  private:
    uint64_t p {0};
    uint64_t w {0};
  };
};


#ifdef WITH_XXHASH
template<>
struct KmerHashers<1>
{
  static std::string name() {return "KmerHashers<1> - XXHASH";}
  template<size_t MAX_K>
  struct Hasher : public IKHasher<MAX_K>
  {
    static std::string name()
    {
      return "KmerHashers<1>::Hasher<MAX_K=" + std::to_string(MAX_K) + ">";
    }
    uint64_t operator()(const Kmer<MAX_K>& kmer, uint64_t seed = 0) const final
    {
      return XXH64(kmer.get_data64(), kmer.m_n_data * sizeof(uint64_t), seed);
    }
  };

  template<size_t MAX_K>
  struct WinHasher : public IKHasher<MAX_K>
  {
    WinHasher(uint64_t p, uint64_t w) : p(p), w(w) {}
    static std::string name()
    {
      return "KmerHashers<1>::WinHasher<MAX_K=" + std::to_string(MAX_K) + ">";
    }

    uint64_t operator()(const Kmer<MAX_K>& kmer, uint64_t seed = 0) const final
    {
      return (XXH64(kmer.get_data64(), kmer.m_n_data * sizeof(uint64_t), seed) % w) + (w * p);
    }

  private:
    uint64_t p {0};
    uint64_t w {0};
  };
};

#endif

};

