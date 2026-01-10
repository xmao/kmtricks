# Code Review: Protein Support for kmtricks

## Summary

This PR adds protein sequence support to kmtricks alongside the existing DNA support. The implementation introduces a new alphabet abstraction layer and a dedicated protein counting pipeline.

## Overall Assessment

**Quality: Good with some concerns**

The implementation follows sound architectural patterns with a clean separation between DNA and protein paths. However, there are several issues that should be addressed.

---

## Positive Aspects

### 1. Clean Alphabet Abstraction (`alphabet.hpp`:26-64)
- Well-designed enum for `Alphabet::DNA` vs `Alphabet::PROTEIN`
- Helper functions for conversion, bit calculations, and slot sizing
- Cleanly encapsulates alphabet-specific logic

### 2. Proper Kmer Class Generalization (`kmer.hpp`)
- `bits_per_symbol()` now alphabet-aware instead of hardcoded 2 bits
- `rev_comp()` and `canonical()` correctly return self for protein (no complement)
- Shift/mask operations use `m_bits_per_symbol` consistently

### 3. IO Format Versioning (`io_common.hpp`:32-33)
- `KM_IO_VERSION_ALPHABET` ensures forward compatibility
- Alphabet metadata in file headers enables detection of mixed inputs

### 4. Input Validation (`kmer_file.hpp`:211-214)
- Proper alphabet mismatch detection across merged inputs
- Clear error message: "Alphabet mismatch across k-mer inputs"

---

## Improvement Plan

### Critical Issues (Must Fix)

#### 1. Off-by-one Error in `u8from16` (`kff_file.hpp`:184-188)
```cpp
void u8from16(uint8_t b[2], uint16_t u16)
{
  b[1] = (uint8_t)(u16>>=8);  // Bug: shifts before first assignment
  b[0] = (uint8_t)(u16>>=8);  // Always 0
}
```
**Fix:** Replace with:
```cpp
void u8from16(uint8_t b[2], uint16_t u16)
{
  b[1] = (uint8_t)u16;
  b[0] = (uint8_t)(u16 >> 8);
}
```

#### 2. Missing Protein Test Cases (`kmer_test.cpp`)
No protein-specific tests were added. The test file only tests DNA sequences.

**Action:** Add tests for:
- Protein kmer encoding/decoding
- 5-bit packing correctness
- `rev_comp()` returning self for protein
- Alphabet mismatch detection

---

### Moderate Issues (Should Fix)

#### 3. Redundant Mask Check (`protein.hpp`:84)
```cpp
mask.rem_mask = (rem_bits == 64) ? ~0ULL : ((1ULL << rem_bits) - 1);
```
If `rem_bits == 64`, then `has_rem` would be false (since `rem_bits = used_bits % 64`), so this branch is unreachable.

**Action:** Remove the unreachable ternary or add a comment explaining the defensive check.

#### 4. Magic Number for Protein Symbol Count (`kmer.hpp`:73)
```cpp
const char bToP[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'};
```
**Action:** Add a named constant `PROTEIN_ALPHABET_SIZE = 21`.

#### 5. Hardcoded Partition Minimum (`protein.hpp`:114-121)
```cpp
inline uint32_t normalize_nb_parts(uint32_t nb_parts)
{
  if (nb_parts == 0)
    nb_parts = 4;
  if (nb_parts < 4)
    nb_parts = 4;
  return nb_parts;
}
```
**Action:** Document why minimum of 4 partitions exists or make it configurable.

#### 6. Missing `clear` Parameter in Protein Count (`protein.hpp`:612-625)
The `ProteinCountTask` doesn't support the `clear` parameter that the DNA path supports.

**Action:** Add `clear` parameter support for feature parity.

---

### Minor Issues (Nice to Fix)

#### 7. Namespace Typo (`kmer.hpp`:813, `kmer_file.hpp`:392, `kff_file.hpp`:411)
```cpp
};  // namespace kmdiff
```
**Fix:** Change to `// namespace km`.

#### 8. Histogram Hardcoded Range (`protein.hpp`:495)
```cpp
hists[i] = opt->hist ? std::make_shared<KHist>(i, opt->kmer_size, 1, 255) : nullptr;
```
**Action:** Define a constant for histogram range (1-255).

#### 9. VLA Usage (`kff_file.hpp`:236, 256, 265)
```cpp
uint8_t minim[minimizer.size()];  // VLA - not standard C++
uint8_t seq[superk.size()];
```
**Action:** Replace VLAs with `std::vector` for C++ standard compliance and portability.

---

### Design Considerations (Future Work)

#### 10. Separate Protein Pipeline vs Unified Path
The current design creates a parallel protein-specific pipeline (`run_protein_all`, `run_protein_count`, etc.). While this works, it creates code duplication in the merge/aggregate paths.

**Future consideration:** Template-based dispatching could reduce duplication.

#### 11. KFF Format for Protein
The implementation creates a custom "KFFP" format (`KFFP_MAGIC = 0x504B4D4B`) rather than extending KFF.

**Action:** Document the decision rationale in user-facing documentation.

#### 12. Superk Skipped for Protein (`protein.hpp`:648-653)
```cpp
template<size_t MAX_K>
inline void run_protein_superk(superk_options_t opt)
{
  spdlog::info("Protein super-k-mer generation is not used; skipping.");
}
```
Protein doesn't use super-k-mers (no minimizers without canonicalization). This is correct behavior.

**Action:** Document this limitation in user-facing docs.

---

## Files Changed Risk Assessment

| File | Changes | Risk |
|------|---------|------|
| `protein.hpp` (new) | New protein pipeline | Medium - core logic |
| `alphabet.hpp` (new) | Alphabet abstraction | Low - simple utilities |
| `kmer.hpp` | Generalized for alphabet | Medium - core data structure |
| `kff_file.hpp` | Added protein KFF format | High - has `u8from16` bug |
| `io_common.hpp` | Added alphabet to headers | Low - additive change |
| `kmer_file.hpp` | Alphabet-aware IO | Low - defensive checks added |
| `cmd.hpp` | Dispatch to protein path | Low - branching logic |
| `cli_common.hpp` | Added `--alphabet` flag | Low - user interface |

---

## Priority Order for Fixes

1. **P0 (Critical):** Fix `u8from16` bug in `kff_file.hpp`
2. **P0 (Critical):** Add protein unit tests
3. **P1 (High):** Replace VLAs with `std::vector`
4. **P1 (High):** Fix namespace typos
5. **P2 (Medium):** Add `clear` parameter to protein count
6. **P2 (Medium):** Define constants for magic numbers
7. **P3 (Low):** Documentation improvements
