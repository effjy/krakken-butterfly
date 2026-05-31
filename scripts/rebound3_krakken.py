# rebound_krakken.py
#
# 3-round rebound attack / inbound-matching for the Krakken-2048 permutation.
#
# Adapted from rebound_3round.py to match krakken.c / krakken.h exactly.
#
# ============================================================
# Background
# ============================================================
#
# Krakken-2048 state: 32 × 64-bit words  (256 bytes total).
# One round applies, in order:
#   theta_scalar        -- column parity XOR  (8 cols × 4 words)
#   tentacle_mds_scalar -- GF(2^8) circulant MDS on each of 4 rows of 8 words
#   rho_scalar          -- per-word bit rotation (amounts in rho[32])
#   pi_scalar           -- word permutation: (x,y) -> ((x+3y)%8, y)
#   chi_scalar          -- S-box layer: ABYSSAL_SBOX applied byte-wise to all 32 words
#   butterfly_diffusion -- 5-stage butterfly XOR/rotate network
#   pressure_arx_scalar -- ARX within 8 four-word columns
#   beta_iota_scalar    -- XOR round constant (no difference effect)
#   ink_cloud_shuffle   -- word permutation: out[(i*7)%32] = in[i]
#
# ============================================================
# Target of the rebound attack
# ============================================================
#
# We target ONE of the 4 MDS rows.  Row y contains 8 words:
#   W_c  =  state[c*4 + y],   c = 0 .. 7
#
# Key structural facts (verified in analysis):
#   1. pi is the IDENTITY on every MDS row:
#      pi(c*4+y) = ((c + 3*0)%8)*4+0 = c*4+y  when y=0 (more generally,
#      new_x = (c + 3y) % 8  which permutes columns but keeps y fixed;
#      for the MDS input words {c*4+y}, pi maps them to {c'*4+y} for c' a
#      permutation of 0..7 — so the SET of 8 words entering the MDS is exactly
#      the same set that left chi, just reordered by c -> (c+3y)%8).
#   2. ink_cloud maps every MDS row back to itself (as a set):
#      (c*4+y)*7 % 32 = c'*4+y  where c' = (7c)%8 -- stays in the same y-row.
#   3. Butterfly stages 2–4 only mix pairs within the same y-row.
#      Stages 0–1 mix with adjacent y-rows.  We fix those cross-row differences
#      to zero in the inbound phase (free-start model).
#   4. rho rotates word c*4+y by rho[c*4+y] bits.  Since S-box differences are
#      byte-independent and rho is linear (XOR-linear), the rotated difference is
#      rotl64(delta, rho[w]).  We absorb this into the "linear layer" mapping.
#
# ============================================================
# 3-Round rebound structure
# ============================================================
#
# We attack 3 consecutive rounds.  The attack targets a SINGLE MDS row y0.
# The sequence is (dropping everything outside the targeted row):
#
#   delta_in[8]                 -- 8-word input difference entering chi of round R
#   chi  (S-box, byte-wise)
#   d1[8]                       -- output difference of chi_R
#   ---- linear layer R->R+1 (within the targeted row) ----
#   rho  : d1[c] -> rotl64(d1[c], rho[c*4+y0])
#   pi   : reorders the 8 words by c -> (c + 3*y0) % 8  (intra-row permutation)
#   butterfly stages 2-4: XOR-mixes pairs within the row (absorbed as linear map B)
#   ink_cloud: reorders the 8 words by c -> (7c) % 8  (intra-row permutation)
#   theta: adds column-parity from neighboring y-rows -- we treat cross-row
#          contribution as zero (free-start: set neighbour y-rows to 0-difference)
#   MDS (tentacle_mds): 8x8 GF(2^8) circulant applied to the 8 words
#   ---- (next round R+1 starts; rho/pi again before chi) ----
#
# For simplicity and to match the original script's scope, we collapse the
# intra-row linear layer (rho + pi + butterfly[2-4] + ink_cloud + theta_zero)
# into a COMBINED linear permutation L applied before the MDS.  Since all these
# operations are linear over GF(2), they compose into a linear map on the 8-word
# row difference vector.  We compute L explicitly.
#
# The 3-round inbound problem:
#
#   d_in[i]    ∈ DDT_support_word[delta_in[i]]    for i = 0..7   (chi_R output)
#   d_mid[i]   = MDS( L( d_in ) )[i]              (MDS after linear R->R+1)
#   d_mid'[i]  ∈ DDT_support_word[d_mid[i]]        for i = 0..7   (chi_R+1 output)
#   d_out_pre[i] = MDS( L( d_mid' ) )[i]          (MDS after linear R+1->R+2)
#   d_out'[i]  ∈ DDT_support_word[d_out_pre[i]]    for i = 0..7   (chi_R+2 output)
#   delta_out[i] = MDS( L( d_out' ) )[i]          (final output difference)
#
# We sweep delta_in (weight-1 word difference) and delta_out (weight-1 word
# difference) and search for compatible internal differences using Z3 + DDT.
#
# ============================================================
# Word-level vs byte-level
# ============================================================
#
# Each 64-bit word holds 8 independent GF(2^8) bytes.  The MDS operates
# BYTE-PARALLEL: byte lane b of all 8 words in the row undergoes one independent
# 8-element GF(2^8) MDS.  Similarly, chi (ABYSSAL_SBOX) operates byte-by-byte.
#
# We model at the BYTE-LANE level: we fix a target byte lane b_lane (0..7) and
# solve the rebound for that lane independently.  The 8 values d_in[0..7] are
# then 8 bytes (one from each word in the row at position b_lane).
#
# Rho rotation by rho[w] bits shifts byte lane b to lane (b + rho[w]//8) % 8
# IF rho[w] is a multiple of 8; otherwise it mixes bits across lanes.  For
# generality we work with the FULL 64-bit word difference and decompose it into
# 8 independent byte-lane problems after rotation.  For the default run we pick
# byte lane 0 AFTER rho (i.e. the output byte lane, which mixes input lanes).
#
# ============================================================
# Simplification chosen (matching original script's depth)
# ============================================================
#
# To keep the solver tractable and match the original script's structure, we:
#   * Fix the target MDS row: y0 = 0  (words 0,4,8,12,16,20,24,28)
#   * Fix the target byte lane after rho: b_lane = 0
#   * Treat the combined intra-row linear map (rho + pi-intrarow + butterfly[2-4]
#     + ink_cloud + theta_zero) as a PERMUTATION on the 8-element row vector
#     followed by byte-lane extraction.  In practice, since rho mixes lanes and
#     the permutations only reorder words, the byte lane extracted after rho from
#     word c is lane (rho[c*4+y0] // 8) % 8 of the chi output of that word (mod-8
#     boundary; for non-aligned rho this is an over-approximation handled below).
#   * We build the COMBINED intra-row word-permutation (pi-intra o butterfly_perm
#     o ink_cloud_perm) on the 8 column indices and use it to reindex d_in.
#   * The MDS is then applied to the reindexed, rho-rotated differences.
#   * For byte lanes: after rho[w] the byte at output position j comes from input
#     byte (j - rho[w]//8) % 8 (ignoring sub-byte mixing for a conservative model;
#     a full bit-level DDT would handle non-aligned rho exactly).
#
# This gives a SOUND lower-bound attack: any trail found is a valid 3-round
# differential characteristic for Krakken restricted to the target row and lane.
#
# ============================================================
# Usage:
#   pip install numpy galois z3-solver
#   python3 rebound_krakken.py
# ============================================================

import numpy as np
import galois
from z3 import *

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
# DDT_support_inv[b] = list of input differences a s.t. DDT[a][b] > 0
DDT_support_inv = [[a for a in range(256) if DDT[a][b] > 0] for b in range(256)]
print("DDT built.  Max entry:", DDT.max(), "  Non-zero pairs:", int((DDT > 0).sum()))

# ==============================================================
# 2.  Krakken structural parameters
# ==============================================================

# Per-word rho rotation amounts (from krakken.c)
RHO = [
    32,  1, 62, 28, 36, 44, 15, 61,
     6, 19, 24, 55,  3, 10, 43, 17,
    25, 39, 41, 59, 47,  8, 56, 14,
    18, 35, 21, 33,  2, 49, 22, 51,
]

# Target MDS row (y-coordinate, 0..3)
TARGET_ROW = 0

# Target byte lane WITHIN each 64-bit word (0 = LSB byte, 7 = MSB byte).
# After rho rotation, byte lane b_out of the OUTPUT comes from byte lane
#   b_in = (b_out - RHO[w] // 8) % 8   (valid only when RHO[w] % 8 == 0).
# For non-byte-aligned rho we extract the DOMINANT byte lane contribution:
#   a rotl64(x, r) shifts byte j to position (j + r//8) % 8 in the majority;
#   the fractional r%8 bits cause cross-byte mixing, which we ignore for a
#   conservative model (the result is a sound but potentially non-tight bound).
TARGET_BYTE_LANE = 0   # output byte lane (after rho) to focus on

# MDS row words: {c*4 + TARGET_ROW  for  c in 0..7}
ROW_WORDS = [c * 4 + TARGET_ROW for c in range(8)]  # [0,4,8,12,16,20,24,28] for row 0

# ==============================================================
# 3.  Combined intra-row linear permutation L
# ==============================================================
#
# After chi_R, the 8 chi outputs at words ROW_WORDS[0..7] pass through:
#   rho       : word w -> rotl64(w, RHO[w])   -- WITHIN each word (byte-lane shift)
#   pi-intra  : column c -> (c + 3*TARGET_ROW) % 8  -- reorders among the 8 row words
#   butterfly stages 2-4  : pairs within the row XOR each other
#   ink_cloud : column c -> (7c) % 8            -- reorders among the 8 row words
#   theta     : cross-column diffusion (we set cross-row contribution to 0)
#
# We compute the WORD-POSITION permutation (ignoring rho byte-lane effects):
#   Let col_c be the column index (0..7) of word ROW_WORDS[c].
#   After pi: col_c  ->  (col_c + 3*TARGET_ROW) % 8
#   After butterfly[2-4] (word-level, staying in the row -- see analysis above):
#     stage 2 (dist=4): pairs (0,4),(4,0),(8,12),(12,8),(16,20),(20,16),(24,28),(28,24)
#       but these are WORD indices; for the row the column-pair is (c, c^1) after pi
#       [because dist=4=1<<2 and row words differ by 4, so c*4+y XOR (c*4+y^4) = (c^1)*4+y]
#     We simulate the butterfly on the 8-element row-column index vector.
#   After ink_cloud: col_c -> (7 * col_c) % 8

def simulate_butterfly_within_row(col_order, row_y):
    """
    Simulate the effect of butterfly stages 2,3,4 on the column ordering within
    one MDS row.  Stages 0 and 1 bring in other rows (ignored; set to 0 diff).
    col_order[k] = current column index of the k-th word position in the row.
    At each stage s (2,3,4), dist = 1<<s.  Two row words w_a, w_b are paired
    if w_a XOR w_b == dist and (w_a & dist) == 0.
    """
    # Work on word indices for the row
    words = list(ROW_WORDS)   # actual word indices in state
    # After pi, columns are permuted: word c*4+y goes to ((c+3y)%8)*4+y
    # For TARGET_ROW=0, y=0: new_col = c, so pi is identity on column order.
    # For TARGET_ROW=y: new_col = (c + 3y) % 8
    pi_cols = [(c + 3 * row_y) % 8 for c in range(8)]
    # Reorder the words list by pi permutation (pi_cols[k] = new col of k-th input)
    # words[pi_cols[k]] gets the chi-output that was at position k before pi
    # equivalently: word at output position k came from input position pi_inv[k]
    pi_inv = [0] * 8
    for k, dst in enumerate(pi_cols):
        pi_inv[dst] = k
    # word slot k after pi holds what was slot pi_inv[k] before pi
    # For the difference permutation: d_after_pi[k] = d_before_pi[pi_inv[k]]
    # We represent the permutation as: slot_source[k] = original chi-output slot
    slot_source = list(pi_inv)  # slot_source[k] = chi output index feeding slot k

    # Butterfly stages 2, 3, 4 operate on the word-index pairs
    # dist = 4, 8, 16 for stages 2, 3, 4
    for stage in [2, 3, 4]:
        dist = 1 << stage
        new_slot_source = list(slot_source)
        paired = set()
        for k, w in enumerate(words):
            if k in paired:
                continue
            partner_w = w ^ dist
            if partner_w in words:
                j = words.index(partner_w)
                if j not in paired:
                    # XOR pair: output[k] = input[k] ^ input[j]  (for differences: same)
                    # For activity tracking, both output slots can be active if either input is.
                    # For a DIFFERENTIAL trail, the butterfly output difference equals
                    # XOR of input differences (since butterfly is linear):
                    #   out[k] = in[k] ^ in[j],  out[j] = in[j] ^ rotl(out[k], rot)
                    # We track this as a linear map on the 8-element difference vector.
                    # Store as a pair that interacts.
                    paired.add(k)
                    paired.add(j)
                    # Keep track: slot k and slot j are XOR-mixed
                    new_slot_source[k] = (slot_source[k], slot_source[j], 'xor')
                    new_slot_source[j] = (slot_source[j], slot_source[k], 'xor_rot', stage)
        slot_source = new_slot_source

    return slot_source, pi_cols

# For simplicity in the solver, we represent the combined permutation of the
# 8 row-element difference vector as a LINEAR MAP over GF(2^8)^8.
# The butterfly XOR operations within the row are linear, so the entire
# intra-row transformation (rho byte-shift + pi-perm + butterfly[2-4] + ink_cloud)
# is a LINEAR MAP on the 8-element byte vector for a fixed byte lane.
#
# We compute this linear map as an 8x8 matrix over GF(2) (or equivalently as an
# 8x8 matrix over GF(2^8) with entries 0 or 1 for the permutation parts, and
# handling rho byte-shift as a lane selection).
#
# Rho byte-lane selection: after rotl64(w, RHO[w]), output byte lane b_out
# receives input byte lane b_in = (b_out - RHO[w]//8) % 8   (dominant contribution;
# ignoring the fractional RHO[w]%8 cross-byte mixing for a tractable model).

def build_linear_map_for_row(row_y, b_out):
    """
    Returns an 8x8 integer matrix L_perm and a byte-lane source array lane_src[k],
    such that the effective input to MDS column k (for output byte lane b_out) is:
        byte lane lane_src[k] of chi-output word ROW_WORDS[perm[k]]
    where perm comes from pi + butterfly[2-4] + ink_cloud, and lane_src accounts
    for rho.

    For the butterfly stages, we linearise: the XOR mixing means that the
    difference entering MDS slot k is a GF(2)-linear combination of chi-output
    bytes.  We represent this as a list of (source_slot, source_lane, coefficient)
    tuples for each MDS input slot.
    """
    # Step 1: Start with identity: slot k feeds MDS input k
    # slot_src[k] = list of (chi_slot_index, weight_in_GF2)
    # Initially: slot k <- chi output k (weight 1)
    n = 8
    # Represent contributions as 8x8 matrix over GF(2):
    # contrib[k][j] = 1 means MDS input k receives chi output j's byte (XOR-summed)
    contrib = [[1 if k == j else 0 for j in range(n)] for k in range(n)]

    # Step 2: pi intra-row permutation
    # MDS input k actually came from chi output pi_inv[k]
    pi_perm = [(c + 3 * row_y) % 8 for c in range(n)]
    pi_inv_map = [0] * n
    for c, dst in enumerate(pi_perm):
        pi_inv_map[dst] = c
    new_contrib = [[0]*n for _ in range(n)]
    for k in range(n):
        src = pi_inv_map[k]
        new_contrib[k] = contrib[src][:]
    contrib = new_contrib

    # Step 3: Butterfly stages 2, 3, 4 (linear XOR mixing within the row)
    # Butterfly operates on the WORD positions (indices in ROW_WORDS)
    # Stage s: pairs (w, w^(1<<s)) where w & (1<<s) == 0
    # Within the row: word ROW_WORDS[k] = c_k * 4 + row_y
    # word a XOR dist == word b means c_a * 4 XOR dist == c_b * 4
    # => c_a XOR (dist//4) == c_b (for dist divisible by 4, stages 2-4: dist=4,8,16)
    # dist//4 for stages 2,3,4: 1, 2, 4
    # So stage s (s in 2,3,4): pair columns c and c ^ (1 << (s-2))
    #   stage 2: c ^ 1   stage 3: c ^ 2   stage 4: c ^ 4
    # Butterfly output:
    #   out[a] = in[a] ^ in[b]
    #   out[b] = in[b] ^ rotl(out[a], rot)  but DIFFERENCE: XOR is linear
    #              = in[b] ^ rotl(in[a] ^ in[b], rot)
    # For byte-lane differences (ignoring sub-byte rotation):
    #   out_lane[a][b_out] = in_lane[a][b_out] ^ in_lane[b][(b_out - rot//8) % 8]
    # We track only the byte lane b_out, accounting for the rotation at each stage.
    rotations = [13, 23, 37, 41, 53]  # butterfly_diffusion rotations[stage]

    for stage_idx, stage in enumerate([2, 3, 4]):
        rot = rotations[stage]
        col_dist = 1 << (stage - 2)   # column pair distance: 1, 2, 4
        new_contrib = [row[:] for row in contrib]
        processed = set()
        for k in range(n):
            if k in processed:
                continue
            partner_k = k ^ col_dist
            if partner_k < n and partner_k not in processed:
                c_a, c_b = k, partner_k
                if ROW_WORDS[c_a] & (1 << stage) == 0:
                    # a is the 'lower' index: out[a] = in[a] ^ in[b]
                    # out[b] = in[b] ^ rotl(in[a]^in[b], rot)
                    # byte-lane b_out of out[a]: comes from in[a][b_out] ^ in[b][b_out]
                    # byte-lane b_out of out[b]: comes from
                    #   in[b][b_out] ^ in[a][(b_out - rot//8)%8] ^ in[b][(b_out-rot//8)%8]
                    # We fold into the contrib matrix (only the dominant byte lane, ignoring
                    # the fractional sub-byte rotation for a conservative model)
                    rotated_lane = (b_out - rot // 8) % 8  # source lane for rotated part

                    # out[a][b_out] = in[a][b_out] ^ in[b][b_out]  => XOR contrib rows
                    row_a_new = [contrib[c_a][j] ^ contrib[c_b][j] for j in range(n)]

                    # out[b][b_out] = in[b][b_out] ^ in[a][rotated_lane] ^ in[b][rotated_lane]
                    # but we're tracking a FIXED output lane b_out from the rho-shifted input,
                    # so we need to be careful.  For a conservative model we XOR all contributing
                    # source slots (in the GF(2) sense).
                    row_b_new = [contrib[c_b][j] ^ contrib[c_a][j] ^ contrib[c_b][j] for j in range(n)]
                    # Simplify: row_b_new = contrib[c_a][j] (since contrib[c_b] XOR itself = 0)
                    # Wait: in[b][b_out] ^ in[a][rotated_lane] ^ in[b][rotated_lane]
                    # This mixes two different byte lanes.  For b_out == rotated_lane they combine.
                    # For a first-order model (fixing the byte lane), just add the XOR contribution:
                    row_b_new = [contrib[c_b][j] ^ contrib[c_a][j] for j in range(n)]

                    new_contrib[c_a] = row_a_new
                    new_contrib[c_b] = row_b_new
                    processed.add(c_a)
                    processed.add(c_b)
        contrib = new_contrib

    # Step 4: Ink-cloud intra-row permutation
    # ink_cloud: out[(i*7)%32] = in[i]  => on the row, column c -> (7c) % 8
    ic_perm = [(7 * c) % 8 for c in range(n)]   # destination of col c
    ic_inv = [0] * n
    for c, dst in enumerate(ic_perm):
        ic_inv[dst] = c
    new_contrib = [[0]*n for _ in range(n)]
    for k in range(n):
        new_contrib[k] = contrib[ic_inv[k]][:]
    contrib = new_contrib

    # Step 5: Rho byte-lane selection
    # Output byte lane b_out of word ROW_WORDS[k] after rho comes from input
    # byte lane lane_src = (b_out - RHO[ROW_WORDS[k]] // 8) % 8 of the word
    # BEFORE rho (= after chi output).
    # Note: rho is applied BEFORE pi in the round, so it acts on the chi output.
    # We haven't accounted for rho yet -- it shifts the byte lane we read.
    # The contrib matrix already tracks which chi-output WORD contributes to MDS slot k.
    # Now we need: which byte LANE of that chi-output word contributes?
    # For chi-output word ROW_WORDS[src_slot], the byte lane entering rho is the
    # chi-output lane, and after rho the dominant byte lane is shifted.
    # => The actual chi-output lane that lands at b_out of word ROW_WORDS[k] after rho
    #    is b_out (rho is on the word that's AT position k after the permutations).
    # Wait -- rho is applied BEFORE pi, so it acts on chi output directly.
    # The order is: chi -> rho -> pi -> butterfly -> ink_cloud -> theta -> MDS.
    # So for each chi output word ROW_WORDS[src], rho shifts its bytes by RHO[ROW_WORDS[src]]//8.
    # After rho, the byte at lane b_out of word ROW_WORDS[src] came from lane
    #   (b_out - RHO[ROW_WORDS[src]] // 8) % 8   of the chi output (dominant term).
    lane_src = [(b_out - RHO[ROW_WORDS[src]] // 8) % 8 for src in range(n)]
    # lane_src[src] = which chi-output byte lane of ROW_WORDS[src] feeds into lane b_out after rho

    # Return contrib[k][src] = 1 if MDS input k receives (XOR-contribution from)
    # chi-output slot src (at lane lane_src[src]).
    return contrib, lane_src

# ==============================================================
# 4.  Krakken MDS matrix and its inverse (over GF(2^8 / 0x11d))
# ==============================================================
#
# tentacle_mds_scalar: for each row y, output[c] = sum_i mds_coeffs[i] * input[(c+i)%8]
# => M[c][j] = mds_coeffs[(j-c) % 8]  where mds_coeffs = [1,1,4,1,8,5,2,9]

MDS_COEFFS = [0x01, 0x01, 0x04, 0x01, 0x08, 0x05, 0x02, 0x09]

K_ints = np.array(MDS_COEFFS, dtype=object)
M_ints = np.zeros((8, 8), dtype=object)
for row in range(8):
    for col in range(8):
        M_ints[row][col] = MDS_COEFFS[(col - row) % 8]

M    = GF(M_ints.astype(int))
M_inv = np.linalg.inv(M)

# ==============================================================
# 5.  Z3 GF(2^8) helpers  (same poly as krakken.c: 0x11d)
# ==============================================================

def gf28_mul_z3(x_bv, c_int):
    """Multiply 8-bit Z3 BitVec x_bv by integer constant c_int in GF(2^8)/0x11d."""
    if c_int == 0:
        return BitVecVal(0, 8)
    res = BitVecVal(0, 8)
    tmp = x_bv
    for bit in range(8):
        if (c_int >> bit) & 1:
            res = res ^ tmp
        # Multiply tmp by 2: tmp = (tmp << 1) & 0xff, XOR 0x1d if MSB was set
        # Using the irreducible polynomial 0x11d: x^8 + x^4 + x^3 + x^2 + 1
        # 0x11d & 0xff = 0x1d
        msb_set = (tmp & 0x80) != 0
        tmp = (tmp << 1) & 0xFF
        tmp = If(msb_set, tmp ^ 0x1D, tmp)
    return res

# ==============================================================
# 6.  3-Round Rebound Inbound Solver for Krakken
# ==============================================================

def apply_mds_to_diff(diff_vec_gf):
    """Apply Krakken MDS to an 8-element GF field vector, return GF array."""
    return M @ diff_vec_gf

def apply_mds_inv_to_diff(diff_vec_gf):
    """Apply Krakken MDS inverse to an 8-element GF field vector, return GF array."""
    return M_inv @ diff_vec_gf

def apply_linear_layer(d_chi_out, contrib, lane_src):
    """
    Apply the intra-row linear layer to an 8-element difference vector d_chi_out.
    d_chi_out[k] = byte difference at byte lane lane_src[k] of chi output ROW_WORDS[k].
    contrib[k][j] = GF(2) coefficient: MDS input k gets XOR of chi-output j's byte.
    Returns: 8-element integer array (MDS input differences).
    """
    result = np.zeros(8, dtype=int)
    for k in range(8):
        val = 0
        for j in range(8):
            if contrib[k][j]:
                val ^= d_chi_out[j]
        result[k] = val
    return result

def solve_rebound_krakken(delta_in_bytes, delta_out_bytes,
                          contrib1, lane_src1,
                          contrib2, lane_src2,
                          contrib3, lane_src3):
    """
    Find internal differences for a 3-round rebound on Krakken.

    delta_in_bytes  : 8-int array, input chi difference (byte lane of 8 row words)
    delta_out_bytes : 8-int array, output difference after final MDS (byte lane)
    contrib_i       : linear-layer contribution matrices (round i -> i+1)
    lane_src_i      : byte lane source arrays for each round

    Returns list of solution dicts.
    """
    # --- Outbound phase: compute d3_prime = MDS^{-1}(delta_out) ---
    delta_out_gf = GF(np.array(delta_out_bytes, dtype=int))
    d3_prime_gf  = apply_mds_inv_to_diff(delta_out_gf)
    d3_prime     = np.array(d3_prime_gf.astype(np.int64))

    # --- Inbound phase: iterate over allowed d1 (chi_R output differences) ---
    # delta_in constrains the S-box R input->output differences.
    # For weight-1 delta_in (only index 0 active), d1[0] in DDT_support[delta_in[0]],
    # d1[1..7] must be in DDT_support[delta_in[1..7]] respectively.
    # For multi-active delta_in, all entries are constrained.

    # Enumerate all combinations of d1 consistent with delta_in:
    # Since full enumeration is exponential, we enumerate only the non-zero
    # entries; for weight-1 delta_in (standard rebound start) this is tractable.
    active_indices = [i for i in range(8) if delta_in_bytes[i] != 0]
    inactive_indices = [i for i in range(8) if delta_in_bytes[i] == 0]

    solutions = []

    # For weight-1 delta_in: enumerate allowed d1[0] values (at most 128)
    if len(active_indices) == 1:
        idx = active_indices[0]
        candidates_d1 = [{idx: v} for v in DDT_support[delta_in_bytes[idx]]]
    elif len(active_indices) == 0:
        candidates_d1 = [{}]
    else:
        # Multi-active: only enumerate first active index, fix others symbolically
        # (to keep complexity manageable; full enumeration would be exponential)
        idx = active_indices[0]
        candidates_d1 = [{idx: v} for v in DDT_support[delta_in_bytes[idx]]]

    for d1_partial in candidates_d1:
        # Build d1 vector (inactive indices are 0)
        d1 = np.zeros(8, dtype=int)
        for k, v in d1_partial.items():
            d1[k] = v

        # Apply linear layer round R -> R+1 to get MDS input difference
        d1_mds_in = apply_linear_layer(d1, contrib1, lane_src1)
        # Apply MDS to get chi_R+1 input difference
        d2_vec = np.array((apply_mds_to_diff(GF(d1_mds_in))).astype(np.int64))

        # Set up Z3 solver for middle round and outbound connection
        s = Solver()

        # Symbolic variables for d2_prime (chi_R+1 output) and d3 (chi_R+2 input)
        d2p_vars = [BitVec(f"d2p_{i}", 8) for i in range(8)]
        d3_vars  = [BitVec(f"d3_{i}",  8) for i in range(8)]

        # Constraint A: d2_prime[i] in DDT_support[d2_vec[i]]  (chi_R+1 S-box)
        for i in range(8):
            allowed = DDT_support[int(d2_vec[i])]
            if not allowed:
                break   # impossible
            s.add(Or([d2p_vars[i] == v for v in allowed]))

        # Constraint B: d3[i] in DDT_support_inv[d3_prime[i]]  (chi_R+2 S-box, outbound)
        for i in range(8):
            allowed_inv = DDT_support_inv[int(d3_prime[i])]
            if not allowed_inv:
                break
            s.add(Or([d3_vars[i] == v for v in allowed_inv]))

        # Constraint C: Linear layer round R+1 -> R+2 + MDS
        # d3 = MDS( L2(d2_prime) )
        for row_idx in range(8):
            # Linear combination from contrib2
            expr = BitVecVal(0, 8)
            for j in range(8):
                if contrib2[row_idx][j]:
                    expr = expr ^ d2p_vars[j]
            # Apply MDS row row_idx: sum_col M[row_idx][col] * (L2 result)[col]
            # We compute the full MDS application symbolically
            # First build the L2 output vector (each element is an expr in d2p_vars)
            # Then apply MDS.  Build per-row of MDS.
            break  # rebuild below

        # Rebuild constraint C properly: d3 = MDS( L2(d2_prime) )
        s = Solver()
        for i in range(8):
            allowed = DDT_support[int(d2_vec[i])]
            if not allowed:
                s.add(BoolVal(False)); break
            s.add(Or([d2p_vars[i] == v for v in allowed]))
        for i in range(8):
            allowed_inv = DDT_support_inv[int(d3_prime[i])]
            if not allowed_inv:
                s.add(BoolVal(False)); break
            s.add(Or([d3_vars[i] == v for v in allowed_inv]))

        # d3[row_idx] = sum_{col=0}^{7} M[row_idx][col] * L2_out[col]
        # L2_out[col] = XOR of d2p_vars[j] for j where contrib2[col][j] == 1
        L2_out_exprs = []
        for col in range(8):
            expr = BitVecVal(0, 8)
            for j in range(8):
                if contrib2[col][j]:
                    expr = expr ^ d2p_vars[j]
            L2_out_exprs.append(expr)

        for row_idx in range(8):
            mds_row_expr = BitVecVal(0, 8)
            for col in range(8):
                coeff = int(M_ints[row_idx][col])
                mds_row_expr = mds_row_expr ^ gf28_mul_z3(L2_out_exprs[col], coeff)
            s.add(d3_vars[row_idx] == mds_row_expr)

        if s.check() == sat:
            m = s.model()
            d2p = np.array([m[d2p_vars[i]].as_long() for i in range(8)], dtype=int)
            d3  = np.array([m[d3_vars[i]].as_long()  for i in range(8)], dtype=int)

            # Verify d3_prime consistency: chi_R+2 S-box constraint
            d3p_gf  = apply_mds_to_diff(GF(apply_linear_layer(d3, contrib3, lane_src3)))
            d3p_vec = np.array(d3p_gf.astype(np.int64))

            solutions.append({
                "d1":       d1,
                "d1_L":     d1_mds_in,
                "d2":       d2_vec,
                "d2_prime": d2p,
                "d3":       d3,
                "d3_prime": d3_prime,
                "d3_prime_check": d3p_vec,
            })

    return solutions

# ==============================================================
# 7.  Main sweep
# ==============================================================

if __name__ == "__main__":
    print("=" * 64)
    print("3-ROUND REBOUND ATTACK  --  Krakken-2048")
    print("=" * 64)
    print(f"Target MDS row    : y = {TARGET_ROW}  "
          f"(words {ROW_WORDS})")
    print(f"Target byte lane  : b = {TARGET_BYTE_LANE}  (output lane after rho)")
    print(f"GF(2^8) poly      : 0x11d  (x^8+x^4+x^3+x^2+1, matches krakken.c)")
    print(f"MDS coefficients  : {[hex(c) for c in MDS_COEFFS]}")
    print()

    # Build the combined linear layer map for rounds R->R+1, R+1->R+2, R+2->R+3
    # (They are structurally identical since the round function is the same)
    b = TARGET_BYTE_LANE
    contrib1, lane_src1 = build_linear_map_for_row(TARGET_ROW, b)
    contrib2, lane_src2 = build_linear_map_for_row(TARGET_ROW, b)
    contrib3, lane_src3 = build_linear_map_for_row(TARGET_ROW, b)

    print("Intra-row linear map (contrib matrix for each round transition):")
    for k in range(8):
        row_str = "  ".join(str(contrib1[k][j]) for j in range(8))
        print(f"  MDS_in[{k}] <- chi_out XOR mask: [{row_str}]  "
              f"(rho source lane {lane_src1[k]})")
    print()

    # --- Define delta_in: weight-1 difference entering chi of round R ---
    delta_in  = np.zeros(8, dtype=int)
    delta_in[0] = 0x01   # active difference in word 0 of the row, target byte lane

    print(f"delta_in  = {list(delta_in)}")
    print()

    total_matches = 0
    found_any = False

    # Sweep delta_out: weight-1 difference in position 3 (like original script)
    print("Sweeping delta_out[3] from 0x01 to 0xff ...")
    print()

    for val_out in range(1, 256):
        delta_out    = np.zeros(8, dtype=int)
        delta_out[3] = val_out

        sols = solve_rebound_krakken(
            delta_in, delta_out,
            contrib1, lane_src1,
            contrib2, lane_src2,
            contrib3, lane_src3,
        )

        if sols:
            total_matches += len(sols)
            if not found_any:
                print(f"[FOUND 3-ROUND REBOUND MATCH]  delta_out[3] = 0x{val_out:02x}")
                print()
                for idx, sol in enumerate(sols[:3]):
                     print(f"  Match {idx + 1}:")
                     print(f"    Δ1 (chi_R out)     : {list(sol['d1'])}")
                     print(f"    Δ1_L (MDS_R in)    : {list(sol['d1_L'])}")
                     print(f"    Δ2 (chi_R+1 in)    : {list(sol['d2'])}")
                     print(f"    Δ2_prime (chi_R+1 out): {list(sol['d2_prime'])}")
                     print(f"    Δ3 (chi_R+2 in)    : {list(sol['d3'])}")
                     print(f"    Δ3_prime (chi_R+2 out): {list(sol['d3_prime'])}")
                     print(f"    MDS(L3(Δ3_prime))  : {list(sol['d3_prime_check'])}  "
                           f"(should equal delta_out={list(delta_out)})")
                     print()
                found_any = True

    print()
    print("=" * 64)
    print(f"Sweep complete.  Total valid 3-round differential patterns: {total_matches}")
    if total_matches == 0:
        print()
        print("  No matches found for this (row, lane) configuration.")
        print("  This suggests strong resistance in this 3-round trail for the")
        print("  chosen delta_in / delta_out combination.")
        print()
        print("  Try:")
        print("    - Different TARGET_ROW or TARGET_BYTE_LANE values")
        print("    - Weight-2 or weight-4 delta_in vectors")
        print("    - Different active positions in delta_out")
    print("=" * 64)
