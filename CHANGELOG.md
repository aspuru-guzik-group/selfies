# Changelog

## v2.1.1 - 14.07.2022
- Fixed index bug in attribution 

## v2.1.0 - 17.05.2022

### Changed:
- Dropped support for Python 3.5-3.6 and will continue to support only current Python versions.

### Added:
- optional attribution to map encoder/decoder output string back to input string (Issue #48, #79)

## v2.0.0 - 21.10.2021

### Changed:
- Improved SMILES parsing (by using adjacencey lists internally), with tighter error handling
  (e.g. issues #62 and #60).
- Faster and improved kekulization algorithm (issue #55 fixed).
- Support for symbols that are constrained to 0 bonds (e.g., `[CH4]`) or >8 bonds
  (users can now specify custon bond constraints with over 8 bonds).
- New `strict=True` flag to `selfies.encoder`, which raises an error if the input
  SMILES violates the current bond constraints. `True` by default, can be `False` for speed-up (if
  SMILES are guaranteed to be correct).
- Added bond constraints for B (max. 3 bonds) to the default and preset constraints.
- Updated the syntax of SELFIES symbols to be cleaner and more readable.
    - Removing `expl` from atomic symbols, e.g., `[C@@Hexpl]` becommes `[C@@H]`
    - Cleaner branch symbols, e.g., `[BranchL_2]` becomes `[=BranchL]`
    - Cleaner ring symbols, e.g., `[Expl=RingL]` becomes `[=RingL]`
    - If you want to use the old symbols, use the `compatible=True` flag to `selfies.decoder`,
      e.g., `sf.decoder('[C][C][Expl=Ring1]',compatible=True)` (not recommended!)
- More logically consistent behaviour of `[Ring]` symbols.
- Standardized SELFIES alphabet, i.e., no two symbols stand for the same atom/ion (issue #58), e.g.,
  `[N+1]` and `[N+]` are equivalent now.
- Indexing symbols are now included in the alphabet returned by `selfies.get_semantic_robust_alphabet`.

### Removed
- Removed `constraints` flag from `selfies.decoder`; please use `selfies.set_semantic_constraints()`
  and pass in `"hypervalent"` or `"octet_rule"` instead.
- Removed `print_error` flag in `selfies.encoder` and `selfies.decoder`,
  which now raise errors `selfies.EncoderError` and `selfies.DecoderError`, respectively.

### Bug Fixes
- Potential chirality inversion of atoms making ring bonds (e.g. ``[C@@H]12CCC2CCNC1``):
  fixed by inverting their chirality in ``selfies.encoder`` such that they are decoded with
  the original chirality preserved.
- Failure to represent mismatching stereochemical specifications at ring bonds
  (e.g. ``F/C=C/C/C=C\C``): fixed by adding new ring symbols (e.g. ``[-/RingL]``, ``[\/RingL]``, etc.).

---

## v1.0.4 - 23.04.2021
### Added:
 * decoder option for relaxed hypervalence rules, `decoder(...,constraints='hypervalent')`
 * decoder option for strict octet rules, `decoder(...,constraints='octet_rule')`
### Bug Fix:
 * Fixed constraint for Phosphorus

 ---

## v1.0.3 - 13.01.2021
### Added:
 * Support for aromatic Si and Al (is not officially supported by Daylight SMILES, but RDKit supports it and examples exist in PubChem).

 ---

## v1.0.2 - 14.10.2020
### Added:
 * Support for aromatic Te and triple bonds.
 * Inbuild SELFIES to 1hot encoding, and 1hot encoding to SELFIES

### Changed:
 * Added default semantic constraints for charged atoms (single positive/negative charge of `[C]`, `[N]`, `[O]`, `[S]`, `[P]`)
 * Raised the bond capacity of `P` to 7 bonds (from 5 bonds).

### Bug Fixes:
 * Fixed bug: `selfies.decoder` did not terminate for malformed SELFIES
   that are missing the closed bracket `']'`.

---

## v1.0.1 - 25.08.2020
### Changed:
 *  Code so that is compatible with python >= 3.5.
 *  More descriptive error messages.

### Bug Fixes:
 *  Minor bug fixes in the encoder for SMILES ending in branches (e.g. `C(Cl)(F)`),
    and SMILES with ring numbers between branches (e.g. `C(Cl)1(Br)CCCC1`)
 *  Minor bug fix with ring ordering in decoder (e.g. `C1CC2CCC12` vs `C1CC2CCC21`).

---

## v1.0.0 - 17.08.2020:
### Added:
 *  Added semantic handling of aromaticity / delocalization (by kekulizing SMILES with aromatic symbols before
    they are translated into SELFIES by `selfies.encoder`).
 *  Added semantic handling of charged species (e.g. `[CH+]1CCC1`).
 *  Added semantic handling of radical species (`[CH]1CCC1`) or any species with explicit hydrogens (e.g. `CC[CH2]`).
 *  Added semantic handling of isotopes (e.g. `[14CH2]=C` or `[235U]`).
 *  Improved semantic handling of explicit atom symbols in square brackets, e.g. Carbene (`[C]=C`).
 *  Improved semantic handling of chirality (e.g. `O=C[Co@@](F)(Cl)(Br)(I)S`).
 *  Improved semantic handling of double-bond configuration (e.g. `F/C=C/C=C/C`).
 *  Added new functions to the library, such as `selfies.len_selfies` and
    `selfies.split_selfies`.
 *  Added advanced-user functions to the library to customize the SELFIES semantic constraints, e.g.
    `selfies.set_semantic_constraints`. Allows to encode for instance diborane, `[BH2]1[H][BH2][H]1`.
 *  Introduced new padding `[nop]` (no operation) symbol.

### Changed:
 *  Optimized the indexing alphabet (it is base-16 now).
 *  Optimized the behaviours of rings and branches to fix an issue with specific non-standard molecules that could not be translated.
 *  Changed behaviour of Ring/Branch, such that states `X9991-X9993` are not necessary anymore.
 *  Significantly improved encoding and decoding algorithms, it is much faster now.

---

## v0.2.4 - 01.10.2019:
### Added:
 *  Function ``get_alphabet()`` which returns a list of 29 selfies symbols
    whose arbitrary combination produce >99.99% valid molecules.

### Bug Fixes:
 *  Fixed bug which happens when three rings start at one node, and two of
    them form a double ring.
 *  Enabled rings with sizes of up to 8000 SELFIES symbols.
 *  Bug fix for tiny ring to RDKit syntax conversion, spanning multiple
    branches.

We thank Kevin Ryan (LeanAndMean@github), Theophile Gaudin and Andrew Brereton
for suggestions and bug reports.

---

## v0.2.2 - 19.09.2019:

### Added:
 *  Enabled ``[C@],[C@H],[C@@],[C@@H],[H]`` to use in a semantic
    constrained way.

We thank Andrew Brereton for suggestions and bug reports.

---

## v0.2.1 - 02.09.2019:

### Added:
 *  Decoder: added optional argument to restrict nitrogen to 3 bonds.
    ``decoder(...,N_restrict=False)`` to allow for more bonds;
    standard: ``N_restrict=True``.
 *  Decoder: added optional argument make ring-function bi-local
    (i.e. confirms bond number at target).
    ``decoder(...,bilocal_ring_function=False)`` to not allow bi-local ring
    function; standard: ``bilocal_ring_function=True``. The bi-local ring
    function will allow validity of >99.99% of random molecules.
 *  Decoder: made double-bond ring RDKit syntax conform.
 *  Decoder: added state X5 and X6 for having five and six bonds free.

### Bug Fixes:
 * Decoder + Encoder: allowing for explicit brackets for organic atoms, for
   instance ``[I]``.
 * Encoder: explicit single/double bond for non-canonical SMILES input
   issue fixed.
 * Decoder: bug fix for ``[Branch*]`` in state X1.

We thank Benjamin Sanchez-Lengeling, Theophile Gaudin and Zhenpeng Yao
for suggestions and bug reports.

---

## v0.1.1 - 04.06.2019:
 * initial release
