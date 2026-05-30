# 🦑 Krakken-2048 Butterfly

## XOR-Rotation Butterfly Diffusion (XRBD) — A 2048-bit Wide-State Cryptographic Permutation

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Pre--release-orange)]()

---

## 📡 Announcement

> **Source code and verification scripts — coming soon.**

The complete reference implementation of the Krakken Butterfly permutation, along with the full cryptanalysis verification suite (MILP models, SMT proofs, rebound attack scripts, and empirical test harnesses), will be published in this repository shortly.

---

## 📝 Overview

**Krakken Butterfly** is an evolution of the Krakken-2048 Abyssal design, introducing a novel **XOR-Rotation Butterfly Diffusion (XRBD)** layer that replaces the original Ink Cloud shuffle. Key innovations include:

- **XRBD Layer**: A 5-stage logarithmic butterfly network achieving full word-level avalanche across 32 words in a single pass using only XOR and rotation operations.
- **Reduced Round Count**: From 10 rounds (original) to 8 rounds, with improved security margins.
- **Proven Full Connectivity**: Theorem-proven complete dependency graph after five stages — every output word depends on every input word.

The design specification has been submitted to the IACR ePrint Archive and is currently under review.

---

## 📂 Repository Contents (Coming Soon)

```
krakken-butterfly/
├── src/                    # Reference implementation (C/AVX2)
│   ├── krakken_butterfly.c
│   ├── krakken_butterfly.h
│   └── avx2/               # AVX2-optimized routines
├── scripts/                # Cryptanalysis verification suite
│   ├── milp/               # MILP active S-box search models
│   ├── smt/                # Z3 bijectivity proofs
│   ├── rebound/            # 3-round and 4-round rebound attacks
│   ├── avalanche/          # Bit diffusion test harness
│   └── rotational/         # Rotational cryptanalysis tests
├── tests/                  # Unit tests and test vectors
├── benchmarks/             # Performance measurement tools
├── docs/                   # Additional documentation
└── README.md
```

---

## 🔬 Cryptographic Properties

| Property | Value |
|----------|-------|
| State size | 2048 bits (32 × 64-bit words) |
| Rounds | 8 |
| S-box | ABYSSAL (8-bit bijection, Δ = 4) |
| MDS branch number | 9 (optimal) |
| XRBD stages | 5 (distances: 1, 2, 4, 8, 16) |
| Rotation constants | (13, 23, 37, 41, 53) |
| Security margin | ≥ 4 rounds vs. rebound attacks |

---

## 🚀 Getting Started (When Available)

```bash
# Clone the repository
git clone https://github.com/effjy/krakken-butterfly.git
cd krakken-butterfly

# Build the reference implementation
make

# Run verification tests
make test

# Run cryptanalysis suite
make cryptanalysis
```

---

## 📄 Specification

The full design specification is documented in `permutation.tex` (pending ePrint approval). Key sections include:

- XOR-Rotation Butterfly Diffusion (XRBD) — novel contribution
- SPN sublayer (Theta → MDS → Chi)
- PRESSURE ARX mixing layer
- Security analysis and experimental validation

---

## 🔗 Related Projects

- [Krakken-2048 Abyssal](https://github.com/effjy/krakken) — Original design (10 rounds)
- [Krakken-Disk](https://github.com/effjy/krakken-disk) — Post-quantum encrypted volume manager
- [Virtual Wipe Turbo](https://github.com/effjy/vwipe) — Forensic-grade data sanitization

---

## 📬 Contact

| Platform | Link |
|----------|------|
| GitHub | [@effjy](https://github.com/effjy) |
| X | [@jfclachance](https://x.com/jfclachance) |
| ORCID | [0009-0005-6377-1675](https://orcid.org/0009-0005-6377-1675) |

---

## 📜 License

MIT License — see [LICENSE](LICENSE) file for details.

---

*🦑 Released into the abyss — 2026*

*"Three paradigms. One round. Total diffusion."*
```
