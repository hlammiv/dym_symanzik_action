# Groups

Multiplication tables loaded at runtime by `load_group()` (see `../group.c`),
plus the generator scripts that produce them (`generators/`).

## File format (what `dym-mod-metro` requires)

Whitespace-separated tokens, in this exact order:

    P                      # group order (1 integer)
    ReTr[0] .. ReTr[P-1]   # P reals   (real part of -Tr of each element)
    ImTr[0] .. ImTr[P-1]   # P reals   (imag part of -Tr of each element)
    mult[0][0] .. mult[P-1][P-1]   # P*P integers, the Cayley table

i.e. `1 + 2P + P²` tokens. `load_group()` reads exactly that many and **ignores
any trailing tokens**, so files with extra data appended (e.g. `mys1080-v4` has
the raw matrices after the table) still load correctly. A file in the *old*
format (order + one Tr row + table = `1 + P + P²` tokens, **no ImTr row**) will
NOT load — it fails the identity/inverse check. The order must satisfy
`P ≤ PMAX` (currently **5040**, set in `../group.h`).

Verify any file with `../verify_group.c` (loads it with the real parser and
checks identity, inverses, in-range entries, and associativity):

    gcc -O2 -std=gnu99 verify_group.c group.o -o verify_group
    ./verify_group groups/six_720x4

## Available groups (all verified)

| File           | Order | Family / gauge group                         |
|----------------|------:|----------------------------------------------|
| `Z2`..`Z16`    |  2–16 | Cyclic Zₙ                                     |
| `Z100`         |   100 | Cyclic                                        |
| `D4`,`D4_old`  |     8 | Dihedral (two bases), 2×2                     |
| `D6`           |    12 | Dihedral                                      |
| `D8`,`D10`     | 16,20 | Dihedral                                      |
| `D16`,`D32`    | 32,64 | Dihedral                                      |
| `Q8`           |     8 | Quaternion                                    |
| `myBT`         |    24 | Binary tetrahedral (SU(2))                    |
| `myBO`         |    48 | Binary octahedral (SU(2))                     |
| `myBI`         |   120 | Binary icosahedral (SU(2))                    |
| `S108`         |   108 | SU(3) subgroup                                |
| `S216`         |   216 | SU(3) subgroup                                |
| `S648`         |   648 | SU(3) subgroup                                |
| `mys1080-v4`   |  1080 | SU(3) subgroup S(1080)                        |
| `sii_60`       |    60 | SU(4) subgroup                                |
| `si_60x4`      |   240 | SU(4) subgroup                                |
| `svii_120x4`   |   480 | SU(4) subgroup                                |
| `sviii_120x4`  |   480 | SU(4) subgroup                                |
| `siii_360x4`   |  1440 | SU(4) subgroup                                |
| `six_720x4`    |  2880 | SU(4) subgroup                                |
| `siv_5040`     |  5040 | SU(4) subgroup *(not tracked in git — 116 MB)*|

## Generators (`generators/`)

Python scripts that build the tables above. SU(2)/SU(3)/cyclic/dihedral
generators (`genZ.py`, `genDn.py`, `genQ8.py`, `genB*.py`, `genS*.py`,
`genClifford*.py`, `genA6.py`) come from the `dym` tree; SU(4) generators
(`gens60.py`, `gens60x4.py`, `gens120x4{1,2}.py`, `gens1440.py`, `gens720x4.py`,
`genstowerpower*.py`, `gens7f.py`) come from the Assi–Lamm subduction toolkit
(`QC/Subduction/SU4/`).

Each generator prints progress lines `Elements: N` while closing the group,
then prints the table in the format above. **To produce a loadable file, strip
those progress lines:**

    python generators/gens720x4.py | grep -v 'Elements:' > groups/six_720x4

(Verified: stripping `Elements:` from the raw `gens720x4` output reproduces the
clean `six_720x4` byte-for-byte.)

### SU(4) provenance notes

- `siv_5040` (order 5040) is the complete table; the table-build loop in the
  generators is O(order³) in Python, so regenerating it is very slow. The clean
  copy here was taken from `lenore:~/Desktop/su4/siv_5040` (complete: 5040 full
  rows) rather than regenerated. It is **git-ignored** because it exceeds
  GitHub's 100 MB file limit — regenerate or re-copy from `lenore` as needed.
- `svi_tower` (order 51840) exists on `lenore` but is truncated (incomplete
  table) and exceeds `PMAX`, so it is not included.
