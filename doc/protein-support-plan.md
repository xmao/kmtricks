# DNA + Protein Support Plan

## Scope
- Support both DNA and protein sequences across all modes: `all`, `repart`, `superk`, `count`, `merge`, `dump`, `aggregate`, `filter`, `combine`.
- Support all data types: `kmer`, `hash`, `vector`, `kff`.
- Maintain backward compatibility for existing DNA runs and files.

## Design Principles
- Default behavior remains DNA unless explicitly set.
- Alphabet is a first-class configuration (runtime flag + persisted in file headers).
- Avoid DNA-specific assumptions in core algorithms; encapsulate alphabet logic.

## Workstreams

### 1) Alphabet & Configuration
- Add an alphabet definition layer in `include/kmtricks/alphabet.hpp` (DNA, PROTEIN).
- Introduce `--alphabet dna|protein` in `src/cli.cpp`.
- Thread alphabet into `include/kmtricks/cmd.hpp`, `include/kmtricks/task_scheduler.hpp`, `include/kmtricks/task.hpp`, and `src/kmtricks.cpp`.

### 2) Core K-mer Representation
- Refactor `include/kmtricks/kmer.hpp` to use alphabet-driven `bits_per_symbol`.
- Update `include/kmtricks/minimizer.hpp`, `include/kmtricks/superk.hpp`, `include/kmtricks/repartition.hpp`, and `include/kmtricks/utils.hpp` to remove hard-coded `2` bit shifts and masks.
- Make `rev_comp`/`canonical` no-op for protein.

### 3) Pipeline Parity (All Modes)
- Keep the existing GATB DNA path under `include/kmtricks/gatb/*`.
- Add a protein backend (new module set) for parsing, counting, superk/repartition, and sorting.
- In `include/kmtricks/task.hpp`, dispatch per mode to DNA or protein implementations.

### 4) IO Formats & Compatibility
- Add alphabet metadata to headers in `include/kmtricks/io/kmer_file.hpp`, `include/kmtricks/io/matrix_file.hpp`, `include/kmtricks/io/pa_matrix_file.hpp`, `include/kmtricks/io/hist_file.hpp`, and `include/kmtricks/io/io_common.hpp`.
- Ensure merge/aggregate paths reject mixed alphabet inputs with clear errors.
- For `include/kmtricks/io/kff_file.hpp`, extend encoding for protein or define a kmtricks-specific variant.

### 5) Hash/Vector/Plugin Parity
- Ensure hashing uses packed bytes independent of alphabet in `include/kmtricks/kmer_hash.hpp`.
- Add alphabet metadata to plugins in `include/kmtricks/plugin.hpp` and initialize in `include/kmtricks/cmd.hpp`.

### 6) Tests & Docs
- Add protein test cases under `tests/` (k-mer encode/decode, IO round-trip, mode smoke tests).
- Update `README.md` and `doc/` for the new flag and limitations.

## Compatibility & Acceptance Criteria
- Old DNA files remain readable; new files record alphabet.
- All modes pass end-to-end for DNA and protein.
- Mixed-alphabet merges are rejected with explicit messages.

## Open Questions
- Protein alphabet policy: strict 20 AA vs allowing `B/Z/J/X/*` and unknown handling.
- KFF interoperability: adopt a standard or define kmtricks-specific encoding.
