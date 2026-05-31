# rebound4_krakken.py
#
# 4-round rebound attack / inbound-matching for the Krakken-2048 permutation.
# Optimized version using global solver and pre-defined DDT relation.

import numpy as np
import galois
from z3 import *
import time

# ==============================================================
# 1.  ABYSSAL S-box and DDT
# ==============================================================

GF = galois.GF(2**8, irreducible_poly=galois.Poly.Int(0x11d))

SBOX = np.array([
    0xA5,0xB6,0xDE,0xF7,0x18,0x37,0x8C,0xC1,0x89,0xDA,0x1E,0x85,0x31,0xF0,0x97,0x77,
    0x41,0x14,0xE8,0xC8,0x8A,0x04,0xB5,0x69,0x1D,0x2B,0x0F,0x2C,0x4E,0x19,0xCC,0x79,
    0xD7,0x4D,0x7D,0x43,0x03,0x3A,0x13,0x92,0x32,0xD9,0x75,0xDF,0xAD,0x81,0xC3,0xF1,
    0xF9,0xA7,0xE2,0x35,0x02,0xDD,0x61,0xA2,0x50,0xE1,0x09,0xC5,0xE3,0x71,0xCB,0x99,
    0x9C,0xB1,0x23,0x86,0x3B,0x93,0x24,0xE9,0xF6,0xB4,0x6A,0x66,0xFE,0x7A,0x3E,0x28,
    0x6E,0xF2,0x9B,0xF8,0x3F,0x2A,0x98,0x10,0xA1,0xFB,0x45,0x36,0x64,0x57,0x8F,0x72,
    0x8B,0x29,0x56,0xFD,0xF4,0xA4,0xED,0xA6,0x76,0xEB,0x6B,0x4A,0xC7,0x5E,0x26,0xD0,
    0x5F,0xCA,0x87,0x52,0x01,0x16,0x67,0xB9,0x74,0x4B,0xCF,0xD2,0x60,0x2F,0x49,0x6F,
    0x39,0x1C,0x5D,0x53,0xE6,0x3C,0xC6,0x7F,0xEA,0xE5,0xBE,0x00,0x65,0x88,0x83,0xE4,
    0x0C,0x38,0x2D,0x80,0xB0,0xAB,0x44,0x84,0x08,0x0D,0xB8,0x51,0x9A,0x2E,0x91,0x68,
    0x40,0x0A,0xFC,0x82,0xBA,0xCE,0x0B,0xFA,0x1A,0x5B,0x62,0x22,0xC9,0x3D,0x8D,0x06,
    0x55,0xD5,0x78,0xAE,0x27,0x9D,0x9E,0xAF,0xB7,0x4F,0xDC,0x9F,0x42,0xA3,0xBC,0x15,
    0xB2,0xDB,0x11,0xA9,0x5C,0xE7,0x7B,0xEF,0xFF,0xC2,0x25,0xEE,0x73,0xF5,0xD6,0x48,
    0x4C,0x21,0x70,0xD1,0x30,0x54,0xA0,0xB3,0x94,0x07,0x58,0xAA,0x96,0x1B,0x1F,0x0E,
    0xD8,0x17,0xE0,0xBB,0x46,0x6C,0xAC,0xA8,0x05,0x7E,0x8E,0x33,0xC4,0xD4,0x59,0xBD,
    0xBF,0xF3,0x20,0x34,0x90,0xCD,0xEC,0x63,0x47,0x95,0x12,0x6D,0xD3,0x5A,0xC0,0x7C,
], dtype=int)

# --- DDT over GF(2^8) ---
print("Building DDT ...")
DDT = np.zeros((256, 256), dtype=np.int32)
for inp in range(256):
    s_inp = SBOX[inp]
    for a in range(256):
        b = s_inp ^ SBOX[inp ^ a]
        DDT[a][b] += 1

# DDT_support[a]  = list of output differences b s.t. DDT[a][b] > 0
DDT_support     = [[b for b in range(256) if DDT[a][b] > 0] for a in range(256)]
print("DDT built.  Max entry:", DDT.max(), "  Non-zero pairs:", int((DDT > 0).sum()))

# ==============================================================
# 2.  Krakken structural parameters
# ==============================================================

RHO = [
    32,  1, 62, 28, 36, 44, 15, 61,
     6, 19, 24, 55,  3, 10, 43, 17,
    25, 39, 41, 59, 47,  8, 56, 14,
    18, 35, 21, 33,  2, 49, 22, 51,
]

TARGET_ROW = 0
TARGET_BYTE_LANE = 0
ROW_WORDS = [c * 4 + TARGET_ROW for c in range(8)]

# ==============================================================
# 3.  Combined intra-row linear permutation L
# ==============================================================

def build_linear_map_for_row(row_y, b_out):
    n = 8
    contrib = [[1 if k == j else 0 for j in range(n)] for k in range(n)]

    pi_perm = [(c + 3 * row_y) % 8 for c in range(n)]
    pi_inv_map = [0] * n
    for c, dst in enumerate(pi_perm):
        pi_inv_map[dst] = c
    new_contrib = [[0]*n for _ in range(n)]
    for k in range(n):
        src = pi_inv_map[k]
        new_contrib[k] = contrib[src][:]
    contrib = new_contrib

    rotations = [13, 23, 37, 41, 53]

    for stage_idx, stage in enumerate([2, 3, 4]):
        rot = rotations[stage]
        col_dist = 1 << (stage - 2)
        new_contrib = [row[:] for row in contrib]
        processed = set()
        for k in range(n):
            if k in processed:
                continue
            partner_k = k ^ col_dist
            if partner_k < n and partner_k not in processed:
                c_a, c_b = k, partner_k
                if ROW_WORDS[c_a] & (1 << stage) == 0:
                    row_a_new = [contrib[c_a][j] ^ contrib[c_b][j] for j in range(n)]
                    row_b_new = [contrib[c_b][j] ^ contrib[c_a][j] for j in range(n)]
                    new_contrib[c_a] = row_a_new
                    new_contrib[c_b] = row_b_new
                    processed.add(c_a)
                    processed.add(c_b)
        contrib = new_contrib

    ic_perm = [(7 * c) % 8 for c in range(n)]
    ic_inv = [0] * n
    for c, dst in enumerate(ic_perm):
        ic_inv[dst] = c
    new_contrib = [[0]*n for _ in range(n)]
    for k in range(n):
        new_contrib[k] = contrib[ic_inv[k]][:]
    contrib = new_contrib

    lane_src = [(b_out - RHO[ROW_WORDS[src]] // 8) % 8 for src in range(n)]
    return contrib, lane_src

# ==============================================================
# 4.  Krakken MDS matrix (over GF(2^8 / 0x11d))
# ==============================================================

MDS_COEFFS = [0x01, 0x01, 0x04, 0x01, 0x08, 0x05, 0x02, 0x09]

K_ints = np.array(MDS_COEFFS, dtype=object)
M_ints = np.zeros((8, 8), dtype=object)
for row in range(8):
    for col in range(8):
        M_ints[row][col] = MDS_COEFFS[(col - row) % 8]

M = GF(M_ints.astype(int))

# ==============================================================
# 5.  Z3 GF(2^8) helpers
# ==============================================================

def gf28_mul_z3(x_bv, c_int):
    if c_int == 0:
        return BitVecVal(0, 8)
    res = BitVecVal(0, 8)
    tmp = x_bv
    for bit in range(8):
        if (c_int >> bit) & 1:
            res = res ^ tmp
        msb_set = (tmp & 0x80) != 0
        tmp = (tmp << 1) & 0xFF
        tmp = If(msb_set, tmp ^ 0x1D, tmp)
    return res

def apply_mds_to_diff(diff_vec_gf):
    return M @ diff_vec_gf

def apply_linear_layer(d_chi_out, contrib, lane_src):
    result = np.zeros(8, dtype=int)
    for k in range(8):
        val = 0
        for j in range(8):
            if contrib[k][j]:
                val ^= d_chi_out[j]
        result[k] = val
    return result

# ==============================================================
# 6.  4-Round Rebound Inbound Solver for Krakken
# ==============================================================

def solve_rebound_krakken_4round(s, ddt_func, delta_in_bytes, delta_out_bytes,
                                 contrib1, lane_src1,
                                 contrib2, lane_src2,
                                 contrib3, lane_src3,
                                 contrib4, lane_src4):
    """
    Find internal differences for a 4-round rebound on Krakken.
    """
    active_indices = [i for i in range(8) if delta_in_bytes[i] != 0]

    # Enumerate concrete d1 values matching delta_in
    if len(active_indices) == 1:
        idx = active_indices[0]
        candidates_d1 = [{idx: v} for v in DDT_support[delta_in_bytes[idx]]]
    elif len(active_indices) == 0:
        candidates_d1 = [{}]
    else:
        idx = active_indices[0]
        candidates_d1 = [{idx: v} for v in DDT_support[delta_in_bytes[idx]]]

    solutions = []

    for d1_partial in candidates_d1:
        d1 = np.zeros(8, dtype=int)
        for k, v in d1_partial.items():
            d1[k] = v

        # Compute concrete d2 = MDS( L1(d1) )
        d1_mds_in = apply_linear_layer(d1, contrib1, lane_src1)
        d2_vec = np.array((apply_mds_to_diff(GF(d1_mds_in))).astype(np.int64))

        s.push()

        # Z3 variables for the internal differences
        d2p_vars = [BitVec(f"d2p_{i}", 8) for i in range(8)]
        d3_vars  = [BitVec(f"d3_{i}",  8) for i in range(8)]
        d3p_vars = [BitVec(f"d3p_{i}", 8) for i in range(8)]
        d4_vars  = [BitVec(f"d4_{i}",  8) for i in range(8)]
        d4p_vars = [BitVec(f"d4p_{i}", 8) for i in range(8)]

        # Constraint 1: d2p_vars ∈ DDT_support[d2_vec]
        for i in range(8):
            allowed = DDT_support[int(d2_vec[i])]
            if not allowed:
                s.add(BoolVal(False)); break
            s.add(Or([d2p_vars[i] == v for v in allowed]))

        # Constraint 2: d3_vars = MDS( L2( d2p_vars ) )
        L2_out = []
        for col in range(8):
            expr = BitVecVal(0, 8)
            for j in range(8):
                if contrib2[col][j]:
                    expr = expr ^ d2p_vars[j]
            L2_out.append(expr)

        for row_idx in range(8):
            mds_row_expr = BitVecVal(0, 8)
            for col in range(8):
                coeff = int(M_ints[row_idx][col])
                mds_row_expr = mds_row_expr ^ gf28_mul_z3(L2_out[col], coeff)
            s.add(d3_vars[row_idx] == mds_row_expr)

        # Constraint 3: d3p_vars ∈ DDT_support[d3_vars]
        for i in range(8):
            s.add(ddt_func(d3_vars[i], d3p_vars[i]) == True)

        # Constraint 4: d4_vars = MDS( L3( d3p_vars ) )
        L3_out = []
        for col in range(8):
            expr = BitVecVal(0, 8)
            for j in range(8):
                if contrib3[col][j]:
                    expr = expr ^ d3p_vars[j]
            L3_out.append(expr)

        for row_idx in range(8):
            mds_row_expr = BitVecVal(0, 8)
            for col in range(8):
                coeff = int(M_ints[row_idx][col])
                mds_row_expr = mds_row_expr ^ gf28_mul_z3(L3_out[col], coeff)
            s.add(d4_vars[row_idx] == mds_row_expr)

        # Constraint 5: d4p_vars ∈ DDT_support[d4_vars]
        for i in range(8):
            s.add(ddt_func(d4_vars[i], d4p_vars[i]) == True)

        # Constraint 6: delta_out = MDS( L4( d4p_vars ) )
        L4_out = []
        for col in range(8):
            expr = BitVecVal(0, 8)
            for j in range(8):
                if contrib4[col][j]:
                    expr = expr ^ d4p_vars[j]
            L4_out.append(expr)

        for row_idx in range(8):
            mds_row_expr = BitVecVal(0, 8)
            for col in range(8):
                coeff = int(M_ints[row_idx][col])
                mds_row_expr = mds_row_expr ^ gf28_mul_z3(L4_out[col], coeff)
            s.add(mds_row_expr == int(delta_out_bytes[row_idx]))

        if s.check() == sat:
            m = s.model()
            d2p = np.array([m[d2p_vars[i]].as_long() for i in range(8)], dtype=int)
            d3  = np.array([m[d3_vars[i]].as_long() for i in range(8)], dtype=int)
            d3p = np.array([m[d3p_vars[i]].as_long() for i in range(8)], dtype=int)
            d4  = np.array([m[d4_vars[i]].as_long() for i in range(8)], dtype=int)
            d4p = np.array([m[d4p_vars[i]].as_long() for i in range(8)], dtype=int)

            solutions.append({
                "d1":       d1,
                "d1_L":     d1_mds_in,
                "d2":       d2_vec,
                "d2_prime": d2p,
                "d3":       d3,
                "d3_prime": d3p,
                "d4":       d4,
                "d4_prime": d4p,
            })

        s.pop()

    return solutions

# ==============================================================
# 7.  Main sweep
# ==============================================================

if __name__ == "__main__":
    print("=" * 64)
    print("4-ROUND REBOUND ATTACK  --  Krakken-2048")
    print("=" * 64)
    print(f"Target MDS row    : y = {TARGET_ROW}  "
          f"(words {ROW_WORDS})")
    print(f"Target byte lane  : b = {TARGET_BYTE_LANE}  (output lane after rho)")
    print(f"GF(2^8) poly      : 0x11d  (x^8+x^4+x^3+x^2+1, matches krakken.c)")
    print(f"MDS coefficients  : {[hex(c) for c in MDS_COEFFS]}")
    print()

    # Pre-build global solver and DDT function relation in Z3
    print("Initializing Z3 global solver and DDT function...")
    t_start = time.time()
    s = Solver()
    ddt_func = Function('DDT', BitVecSort(8), BitVecSort(8), BoolSort())
    for a in range(256):
        for b in range(256):
            val = bool(DDT[a][b] > 0)
            s.add(ddt_func(a, b) == val)
    print(f"Z3 initialization took {time.time() - t_start:.2f} seconds.")
    print()

    b = TARGET_BYTE_LANE
    contrib1, lane_src1 = build_linear_map_for_row(TARGET_ROW, b)
    contrib2, lane_src2 = build_linear_map_for_row(TARGET_ROW, b)
    contrib3, lane_src3 = build_linear_map_for_row(TARGET_ROW, b)
    contrib4, lane_src4 = build_linear_map_for_row(TARGET_ROW, b)

    delta_in  = np.zeros(8, dtype=int)
    delta_in[0] = 0x01   # active difference in word 0 of the row, target byte lane

    print(f"delta_in  = {list(delta_in)}")
    print()

    total_matches = 0
    found_any = False

    # Sweep delta_out: weight-1 difference in position 3
    print("Sweeping delta_out[3] from 0x01 to 0xff ...")
    print()

    for val_out in range(1, 256):
        delta_out    = np.zeros(8, dtype=int)
        delta_out[3] = val_out

        sols = solve_rebound_krakken_4round(
            s, ddt_func, delta_in, delta_out,
            contrib1, lane_src1,
            contrib2, lane_src2,
            contrib3, lane_src3,
            contrib4, lane_src4,
        )

        if sols:
            total_matches += len(sols)
            if not found_any:
                print(f"[FOUND 4-ROUND REBOUND MATCH]  delta_out[3] = 0x{val_out:02x}")
                print()
                for idx, sol in enumerate(sols[:3]):
                     print(f"  Match {idx + 1}:")
                     print(f"    Δ1 (chi_R out)        : {list(sol['d1'])}")
                     print(f"    Δ1_L (MDS_R in)       : {list(sol['d1_L'])}")
                     print(f"    Δ2 (chi_R+1 in)       : {list(sol['d2'])}")
                     print(f"    Δ2_prime (chi_R+1 out): {list(sol['d2_prime'])}")
                     print(f"    Δ3 (chi_R+2 in)       : {list(sol['d3'])}")
                     print(f"    Δ3_prime (chi_R+2 out): {list(sol['d3_prime'])}")
                     print(f"    Δ4 (chi_R+3 in)       : {list(sol['d4'])}")
                     print(f"    Δ4_prime (chi_R+3 out): {list(sol['d4_prime'])}")
                     print()
                found_any = True

    print()
    print("=" * 64)
    print(f"Sweep complete.  Total valid 4-round differential patterns: {total_matches}")
    print("=" * 64)
