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

#include <kmtricks/exceptions.hpp>

namespace km {

enum class Alphabet : uint8_t
{
  DNA = 0,
  PROTEIN = 1
};

inline Alphabet alphabet_from_string(const std::string& value)
{
  if (value == "dna")
    return Alphabet::DNA;
  if (value == "protein")
    return Alphabet::PROTEIN;
  throw ConfigError("Unknown alphabet: " + value);
}

inline std::string alphabet_to_string(Alphabet alphabet)
{
  return alphabet == Alphabet::DNA ? "dna" : "protein";
}

inline uint8_t alphabet_bits(Alphabet alphabet)
{
  return alphabet == Alphabet::DNA ? 2 : 5;
}

inline Alphabet alphabet_from_id(uint8_t id)
{
  return id == 0 ? Alphabet::DNA : Alphabet::PROTEIN;
}

inline uint32_t kmer_slots(uint32_t kmer_size, Alphabet alphabet)
{
  const uint32_t bits = kmer_size * alphabet_bits(alphabet);
  return static_cast<uint32_t>((bits + 63) / 64);
}

}  // namespace km
