# Changelog

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
