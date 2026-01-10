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
#include <cstdint>
#include <string>
#include <memory>
#include <optional>
#include <fstream>
#include <vector>

// ext
#include <kff_io.hpp>

// int
#include <kmtricks/kmer.hpp>
#include <kmtricks/alphabet.hpp>

namespace km {

using kff_t = std::unique_ptr<Kff_file>;
using kff_raw_t = std::unique_ptr<Section_Raw>;
using kff_min_t = std::unique_ptr<Section_Minimizer>;

constexpr uint32_t KFFP_MAGIC = 0x504B4D4B; // "KMKP"
constexpr uint32_t KFFP_VERSION = 1;

template<size_t MAX_C>
class KffWriter
{
public:
  KffWriter(const std::string& path, size_t kmer_size, Alphabet alphabet = Alphabet::DNA)
    : m_kmer_size(kmer_size), m_alphabet(alphabet)
  {
    if (alphabet == Alphabet::PROTEIN)
    {
      m_is_protein = true;
      m_bits_per_symbol = alphabet_bits(alphabet);
      m_bytes_per_kmer = (m_kmer_size * m_bits_per_symbol + 7) / 8;
      m_protein_out.open(path, std::ios::binary | std::ios::out);
      if (!m_protein_out.good())
        throw std::runtime_error("Unable to open " + path);
      m_protein_out.write(reinterpret_cast<const char*>(&KFFP_MAGIC), sizeof(KFFP_MAGIC));
      m_protein_out.write(reinterpret_cast<const char*>(&KFFP_VERSION), sizeof(KFFP_VERSION));
      m_protein_out.write(reinterpret_cast<const char*>(&m_kmer_size), sizeof(m_kmer_size));
      uint32_t count_size = sizeof(typename selectC<MAX_C>::type);
      m_protein_out.write(reinterpret_cast<const char*>(&count_size), sizeof(count_size));
      uint8_t alphabet_id = static_cast<uint8_t>(alphabet);
      m_protein_out.write(reinterpret_cast<const char*>(&alphabet_id), sizeof(alphabet_id));
      m_protein_out.write(reinterpret_cast<const char*>(&m_bits_per_symbol), sizeof(m_bits_per_symbol));
    }
    else
    {
      m_kff_file = std::make_unique<Kff_file>(path, "w");
      uint8_t encoding[] = {0, 1, 3, 2};
      m_kff_file->write_encoding(encoding);

      Section_GV sgv(m_kff_file.get());
      sgv.write_var("k", m_kmer_size);
      sgv.write_var("kmtricks_alphabet", static_cast<uint64_t>(alphabet));
      sgv.write_var("kmtricks_alphabet_bits", static_cast<uint64_t>(alphabet_bits(alphabet)));
      sgv.write_var("max", 1);
      sgv.write_var("data_size", sizeof(typename selectC<MAX_C>::type));
      sgv.close();
      m_kff_sec = std::make_unique<Section_Raw>(m_kff_file.get());
    }
  }

  template<size_t MAX_K>
  void write(const Kmer<MAX_K>& kmer, typename selectC<MAX_C>::type count)
  {
    uint8_t encoded[1024];
    uint8_t counts[sizeof(count)];
    if constexpr(sizeof(count) == 1)
      counts[0] = count;
    else if constexpr(sizeof(count) == 2)
      u8from16(counts, count);
    else
      u8from32(counts, count);

    m_kmer = kmer.to_string();
    if (m_is_protein)
    {
      auto packed = encode_protein_sequence(m_kmer);
      m_protein_out.write(reinterpret_cast<const char*>(packed.data()), packed.size());
      m_protein_out.write(reinterpret_cast<const char*>(counts), sizeof(count));
    }
    else
    {
      encode_sequence(encoded);
      m_kff_sec->write_compacted_sequence(encoded, m_kmer_size, counts);
    }
  }

  void close()
  {
    if (m_is_protein)
    {
      m_protein_out.close();
    }
    else
    {
      m_kff_sec->close();
      m_kff_file->close();
    }
  }

  uint8_t uint8_packing(std::string sequence)
  {
    size_t size = sequence.size();
    uint8_t val = 0;
    for (size_t i = 0; i < size; i++)
    {
      val <<= 2;
      val += (sequence[i] >> 1) & 0b11;
    }
    return val;
  }

  void encode_sequence(uint8_t* encoded)
  {
    size_t size = m_kmer.length();
    size_t remnant = size % 4;
    if (remnant > 0)
    {
      encoded[0] = uint8_packing(m_kmer.substr(0, remnant));
      encoded += 1;
    }

    size_t nb_uint_needed = size / 4;
    for (size_t i = 0; i < nb_uint_needed; i++)
    {
      encoded[i] = uint8_packing(m_kmer.substr(remnant + 4 * i, 4));
    }
  }

private:
  std::vector<uint8_t> encode_protein_sequence(const std::string& sequence) const
  {
    std::vector<uint8_t> packed(m_bytes_per_kmer, 0);
    uint64_t buffer = 0;
    uint8_t bits = 0;
    size_t out_idx = 0;
    for (char c : sequence)
    {
      buffer = (buffer << m_bits_per_symbol) | protein_to_code(static_cast<unsigned char>(c));
      bits += m_bits_per_symbol;
      while (bits >= 8 && out_idx < packed.size())
      {
        bits -= 8;
        packed[out_idx++] = static_cast<uint8_t>((buffer >> bits) & 0xFF);
      }
    }
    if (bits && out_idx < packed.size())
    {
      packed[out_idx] = static_cast<uint8_t>((buffer << (8 - bits)) & 0xFF);
    }
    return packed;
  }

  void u8from32(uint8_t b[4], uint32_t u32)
  {
    b[3] = (uint8_t)u32;
    b[2] = (uint8_t)(u32>>=8);
    b[1] = (uint8_t)(u32>>=8);
    b[0] = (uint8_t)(u32>>=8);
  }

  void u8from16(uint8_t b[2], uint16_t u16)
  {
    b[1] = (uint8_t)u16;
    b[0] = (uint8_t)(u16 >> 8);
  }

private:
  kff_t m_kff_file {nullptr};
  kff_raw_t m_kff_sec {nullptr};
  size_t m_kmer_size;
  std::string m_kmer;
  bool m_is_protein {false};
  Alphabet m_alphabet {Alphabet::DNA};
  uint8_t m_bits_per_symbol {2};
  size_t m_bytes_per_kmer {0};
  std::ofstream m_protein_out;
};

template<size_t MAX_C>
using kff_w_t = std::shared_ptr<KffWriter<MAX_C>>;

using kff_reader_t = std::unique_ptr<Kff_reader>;

template<size_t MAX_C>
class KffSkWriter
{
public:
  KffSkWriter(const std::string& path, size_t kmer_size, size_t minim_size,
              Alphabet alphabet = Alphabet::DNA)
    : m_kmer_size(kmer_size), m_minim_size(minim_size)
  {
    if (alphabet == Alphabet::PROTEIN)
      throw ConfigError("Protein kff superk output is not supported.");
    m_kff_file = std::make_unique<Kff_file>(path, "w");
    uint8_t encoding[] = {0, 1, 3, 2};
    m_kff_file->write_encoding(encoding);

    Section_GV sgv(m_kff_file.get());
    sgv.write_var("k", m_kmer_size);
    sgv.write_var("m", m_minim_size);
    sgv.write_var("kmtricks_alphabet", static_cast<uint64_t>(alphabet));
    sgv.write_var("kmtricks_alphabet_bits", static_cast<uint64_t>(alphabet_bits(alphabet)));
    sgv.write_var("max", 255);
    sgv.write_var("data_size", 1);
    sgv.close();
  }

  ~KffSkWriter() { if (m_current_section) m_current_section->close(); }
  void new_section(const std::string& minimizer)
  {
    if (m_current_section) m_current_section->close();
    m_current_section = std::make_unique<Section_Minimizer>(m_kff_file.get());
    std::vector<uint8_t> minim(minimizer.size());
    encode_sequence(minimizer, minim.data());
    m_current_section->write_minimizer(minim.data());
    m_current_minim = minimizer;
  }

  uint8_t uint8_packing(std::string sequence)
  {
    size_t size = sequence.size();
    uint8_t val = 0;
    for (size_t i = 0; i < size; i++)
    {
      val <<= 2;
      val += (sequence[i] >> 1) & 0b11;
    }
    return val;
  }

  void write(const std::string& superk, size_t minim_pos, std::vector<uint8_t>& vcount)
  {
    std::vector<uint8_t> seq(superk.size());
    encode_sequence(superk, seq.data());
    m_current_section->write_compacted_sequence(seq.data(), superk.size(), minim_pos, vcount.data());
  }

  void encode_sequence(const std::string& superk, uint8_t* encoded)
  {
    size_t size = superk.length();
    size_t remnant = size % 4;
    if (remnant > 0)
    {
      encoded[0] = uint8_packing(superk.substr(0, remnant));
      encoded += 1;
    }

    size_t nb_uint_needed = size / 4;
    for (size_t i = 0; i < nb_uint_needed; i++)
    {
      encoded[i] = uint8_packing(superk.substr(remnant + 4 * i, 4));
    }
  }

private:
  size_t m_kmer_size;
  size_t m_minim_size;
  kff_t m_kff_file;
  std::string m_current_minim;
  kff_min_t m_current_section {nullptr};
};


template<size_t MAX_C>
using kffsk_w_t = std::shared_ptr<KffSkWriter<MAX_C>>;

// from kff-tools encoding
class KffReader
{
  std::string nt[4] = {"A", "C", "T", "G"};
public:
  KffReader(const std::string& path, size_t kmer_size)
  {
    std::ifstream probe(path, std::ios::binary);
    uint32_t magic = 0;
    if (probe.good())
      probe.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    if (magic == KFFP_MAGIC)
    {
      m_is_protein = true;
      uint32_t version = 0;
      probe.read(reinterpret_cast<char*>(&version), sizeof(version));
      probe.read(reinterpret_cast<char*>(&m_kmer_size), sizeof(m_kmer_size));
      probe.read(reinterpret_cast<char*>(&m_count_size), sizeof(m_count_size));
      uint8_t alphabet_id = 0;
      probe.read(reinterpret_cast<char*>(&alphabet_id), sizeof(alphabet_id));
      m_alphabet = alphabet_from_id(alphabet_id);
      probe.read(reinterpret_cast<char*>(&m_bits_per_symbol), sizeof(m_bits_per_symbol));
      m_bytes_per_kmer = (m_kmer_size * m_bits_per_symbol + 7) / 8;
      m_protein_in = std::move(probe);
      m_protein_buffer.resize(m_bytes_per_kmer);
    }
    else
    {
      m_kff_reader = std::make_unique<Kff_reader>(path);
      m_kmer_size = kmer_size;
      m_data_size = 0;
      for (int i=0; i<256; i++)
      {
        for (int j=0; j<4; j++)
        {
          uint8_t n = (i>>(2*j)) & 0b11;
          m_lookup[i] = nt[n] + m_lookup[i];
        }
      }
    }
  }

  template<size_t MAX_K>
  std::optional<Kmer<MAX_K>> read()
  {
    static Kmer<MAX_K> kmer;
    kmer.set_k(m_kmer_size);
    if (m_is_protein)
    {
      Kmer<MAX_K>::set_alphabet(m_alphabet);
      if (!m_protein_in.read(reinterpret_cast<char*>(m_protein_buffer.data()), m_bytes_per_kmer))
        return std::nullopt;
      if (!m_protein_in.read(reinterpret_cast<char*>(m_count_buffer.data()), m_count_size))
        return std::nullopt;
      return Kmer<MAX_K>(decode_protein());
    }
    else
    {
      if (m_kff_reader->has_next())
      {
        m_kff_reader->next_kmer(m_buffer, m_data);
        return Kmer<MAX_K>(to_string());
      }
      return std::nullopt;
    }
  }

private:
  std::string decode_protein() const
  {
    std::string seq;
    seq.reserve(m_kmer_size);
    uint64_t buffer = 0;
    uint8_t bits = 0;
    for (auto byte : m_protein_buffer)
    {
      buffer = (buffer << 8) | byte;
      bits += 8;
      while (bits >= m_bits_per_symbol && seq.size() < m_kmer_size)
      {
        bits -= m_bits_per_symbol;
        uint8_t code = static_cast<uint8_t>((buffer >> bits) & ((1ULL << m_bits_per_symbol) - 1));
        seq.push_back(code_to_char(Alphabet::PROTEIN, code));
      }
    }
    return seq;
  }

  std::string to_string()
  {
    size_t size = m_kmer_size % 4 == 0 ? m_kmer_size / 4 : m_kmer_size / 4 + 1;
    std::string str_kmer = m_lookup[m_buffer[0]];
    if (m_kmer_size % 4 != 0)
    {
      int trunc = 4 - (m_kmer_size % 4);
      str_kmer = str_kmer.substr(trunc);
    }
    for (size_t i=1; i<size; i++)
    {
      uint8_t v = m_buffer[i];
      str_kmer += m_lookup[v];
    }
    return str_kmer;
  }

  kff_reader_t m_kff_reader {nullptr};
  size_t m_kmer_size{0};
  size_t m_data_size{0};
  uint8_t* m_buffer {nullptr};
  uint8_t* m_data {nullptr};
  std::string m_lookup[256];
  bool m_is_protein {false};
  std::ifstream m_protein_in;
  size_t m_bytes_per_kmer {0};
  uint8_t m_bits_per_symbol {0};
  uint32_t m_count_size {0};
  Alphabet m_alphabet {Alphabet::DNA};
  std::vector<uint8_t> m_protein_buffer;
  std::array<uint8_t, 4> m_count_buffer {};
};

using kff_r_t = std::unique_ptr<KffReader>;
}; // end of namespace km
