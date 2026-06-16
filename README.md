<div align="center">

<a href="https://github.com/effjy/krakken-butterfly/"><img src="titles/krakken-2048-butterfly-title.svg" height="52" alt="Krakken-2048 Butterfly"></a>

### XOR-Rotation Butterfly Diffusion (XRBD) — A 2048-bit Wide-State Cryptographic Permutation

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Active-brightgreen)]()
[![DOI](https://img.shields.io/badge/DOI-10.6084/m9.figshare.32527287-blue)](https://doi.org/10.6084/m9.figshare.32527287)

</div>

---

## 📡 Announcement

> **✅ Source code and verification scripts are now available!**

The complete reference implementation of the Krakken Butterfly permutation, along with the full cryptanalysis verification suite (MILP models, SMT proofs, rebound attack scripts, and empirical test harnesses), is now published in this repository.

---

## 📝 Overview

**Krakken Butterfly** is an evolution of the Krakken-2048 Abyssal design, introducing a novel **XOR-Rotation Butterfly Diffusion (XRBD)** layer. Key innovations include:

- **XRBD Layer**: A 5-stage logarithmic butterfly network achieving full word-level avalanche across 32 words in a single pass using only XOR and rotation operations.
- **Reduced Round Count**: From 10 rounds (original) to 8 rounds, with improved security margins.
- **Proven Full Connectivity**: Theorem-proven complete dependency graph after five stages — every output word depends on every input word.

The design specification has been published as a preprint (see below).

---

## 📄 Paper

**XRBD: A Bijective Butterfly Diffusion Layer with Optimal Dependency Expansion for Wide-State Cryptographic Permutations**  
Jean-François Lachance-Caumartin

[![DOI](https://img.shields.io/badge/DOI-10.6084/m9.figshare.32527287-blue)](https://doi.org/10.6084/m9.figshare.32527287)

📄 **Download the full paper:** [`paper.pdf`](paper.pdf)

The paper presents the XRBD primitive, proves its optimal dependency expansion and bijectivity, and integrates it into the Krakken-2048 Butterfly permutation. It includes formal security proofs, cryptanalytic evaluation (MILP, SMT, rebound attacks), and AVX2 performance benchmarks.

---

## 📂 Repository Structure

| File | Description |
|------|-------------|
| `krakken.h` | Header file with function prototypes, S‑box tables, and round constants. |
| `krakken.c` | Core permutation implementation (scalar, constant‑time). |
| `krakken_multi.c` | AVX2‑vectorized and multi‑threaded version. |
| `verify_krakken.c` | Verification test suite (avalanche, rebound, active S‑box, etc.). |
| `Makefile` | Build configuration (GCC, O3 optimizations, pthread). |
| `scripts/` | Cryptanalysis and formal verification scripts (see below). |

### 🔬 Cryptanalysis Scripts (`scripts/`)

| Script | Purpose |
|--------|---------|
| `milp_active_sboxes.py` | MILP model for active S‑box lower bounds (2 rounds). |
| `milp_krakken.py` | Full MILP differential trail search over multiple rounds. |
| `rebound3_krakken.py` | 3‑round rebound attack search using Z3 SMT solver. |
| `rebound4_krakken.py` | 4‑round rebound attack search (demonstrates security margin). |
| `sbox_algebraic_relations.py` | Computes quadratic algebraic relations of the ABYSSAL S‑box. |
| `z3_bijectivity_proof.py` | Formal Z3 proof of bijectivity for PRESSURE ARX and Chi layers. |

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

## 🚀 Getting Started

```bash
# Clone the repository
git clone https://github.com/effjy/krakken-butterfly.git
cd krakken-butterfly

# Build all variants
make

# Run verification tests and benchmarks
make run

# (Optional) Run cryptanalysis scripts individually
python3 scripts/milp_active_sboxes.py
python3 scripts/rebound3_krakken.py
```

### Build Targets

| Target | Description |
|--------|-------------|
| `krakken_scalar` | Scalar (constant‑time) version |
| `krakken_multi`  | AVX2 + multi‑threaded version |
| `verify_krakken` | Verification and cryptanalysis suite |

---

## 📄 Specification

The full design specification is documented in [`paper.pdf`](paper.pdf) (also available via [Figshare DOI](https://doi.org/10.6084/m9.figshare.32527287)). Key sections include:

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
