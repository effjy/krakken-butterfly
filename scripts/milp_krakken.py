# milp_krakken.py
#
# Formal MILP lower-bound search for active S-boxes in the Krakken permutation.
#
# Krakken is a 2048-bit (32 × 64-bit word) permutation with 10 rounds.
# Each round applies, in order:
#   1. theta_scalar          -- column parity XOR diffusion (8 cols × 4 words)
#   2. tentacle_mds_scalar   -- GF(2^8) circulant MDS on each of 4 rows of 8 words
#   3. rho_scalar            -- per-word bit rotation (no word-level activity change)
#   4. pi_scalar             -- word permutation: (x,y) → ((x+3y)%8, y)
#   5. chi_scalar            -- S-box layer: ABYSSAL_SBOX applied byte-wise to all
#                               32 words (= 256 independent 8-bit S-boxes per round)
#   6. butterfly_diffusion   -- 5-stage butterfly XOR/rotate mixing all 32 words
#   7. pressure_arx_scalar   -- ARX within 8 four-word columns
#   8. beta_iota_scalar      -- XOR round constant (no difference effect)
#   9. ink_cloud_shuffle     -- word permutation: out[(i*7)%32] = in[i]
#
# S-box identification
# --------------------
# The S-box layer (chi) applies custom_sbox8_64 to every word, which is just the
# ABYSSAL_SBOX lookup applied independently to each of the 8 bytes in a 64-bit
# word.  There are 32 × 8 = 256 S-boxes per round.
#
# MILP model
# ----------
# We model activity at the **word level** (32 words).  A word is "active" if at
# least one of its 8 bytes is active.  All sub-word operations (rho rotations,
# 11-bit rotation in ink_cloud, ARX) mix bits within a word, so they cannot
# decrease word-level activity; we track them conservatively.
#
# Between consecutive chi applications the linear layer is:
#
#   chi_out  →  butterfly  →  ARX  →  ink_cloud  →  theta  →  MDS  →  rho  →  pi  →  chi_in
#
# We introduce auxiliary word-activity variables at each intermediate stage:
#
#   x[r][w]  : word w active entering chi (S-box input) of round r
#   a[r][w]  : word w active exiting  chi (S-box output)  of round r
#   bf[r][s][w]: word w active after butterfly stage s      of round r
#   t[r][w]  : word w active after ink_cloud               of round r
#   v[r][w]  : word w active after theta                   of round r
#   (MDS constraints on rows of v and x of next round)
#
# Key constraints
# ---------------
# 1. S-box (chi): bijection at word level → a[r][w] == x[r][w].
#    We model this as a[r][w] >= x[r][w] (only direction needed for lower bound;
#    in a bijection output activity = input activity, but for the MILP objective
#    we only count x[r][w], so we just need the propagation direction).
#
# 2. Butterfly (5 stages of XOR-pair diffusion):
#    At each stage s, pairs (i, i^(1<<s)) for (i & (1<<s)) == 0 are mixed.
#    If either input word of a pair is active, the output CAN be active for both.
#    For a lower-bound model: output[w] >= input[w] and output[w] >= input[partner].
#    This is conservative (allows the trail to not spread if there's no reason).
#
# 3. ARX: mixes within 4-word columns (same grouping as theta columns).
#    Conservative: ARX output word w active if any word in the same ARX column
#    was active entering ARX.  We share the butterfly output → ink_cloud transition.
#    (ARX and ink_cloud are merged: we just use bf[r][4][w] as the pre-ink-cloud
#    state since both ARX and ink_cloud only spread activity.)
#
# 4. Ink-cloud shuffle: pure word permutation out[(i*7)%32] = in[i].
#    t[r][(w*7)%32] = bf[r][4][w]  (exact, no new activity).
#
# 5. Theta: column parity XOR diffusion.
#    Columns: col c = {4c, 4c+1, 4c+2, 4c+3}.
#    d[c] = rotr(parity[c-1], 1) XOR parity[c+1].
#    Word 4c+y ^= d[c], so column c is activated if col c-1 or col c+1 is active.
#    For a MILP lower bound: v[r][w] >= t[r][w], and
#    if any word in col c-1 or c+1 is active then all words in col c may become active.
#    We use indicator uc[r][c] = 1 if any word in col c is active (after ink_cloud),
#    then v[r][4c+y] >= uc[r][(c-1)%8] (spread from left neighbor)
#         v[r][4c+y] >= uc[r][(c+1)%8] (spread from right neighbor).
#    This gives a valid lower bound (we are only ADDING activity, never removing).
#
# 6. MDS (tentacle_mds): applied to each of the 4 rows of 8 words.
#    MDS coefficients: [1,1,4,1,8,5,2,9] over GF(2^8) (8×8 circulant).
#    Branch number over GF(2^8): 9  (MDS property; 8×8 MDS → BN = 8+1 = 9).
#    The MDS operates byte-by-byte across the 8 words of a row.
#    Word-level branch number: also 9  (if any byte of a word is active, the word
#    is active; the MDS constraint then applies at word level with the same BN).
#    Constraint: for row y (words {0*4+y, ..., 7*4+y}):
#      let  in_words  = v[r][c*4+y]  for c in 0..7   (MDS input)
#      let  out_words = x[r+1][pi_of(c*4+y)] for c in 0..7  (MDS output → chi input)
#      sum(in_words) + sum(out_words) >= 9 * mds_act[r][y]
#      in_words[c]  >= mds_act[r][y]   -- if row active, each input may be active
#      out_words[c] >= mds_act[r][y]   -- (these are redundant but tighten the LP)
#    Note: "out_words" are the words that the MDS output feeds into, after rho (no
#    position change) and pi (word permutation), entering chi of round r+1.
#
# 7. Pi: word permutation  out[new_x*4+y] = in[x*4+y]  where new_x = (x+3y)%8.
#    This is applied after MDS output → rho → pi → chi_in of next round.
#    We precompute: for each MDS-output position w, its chi_in position.
#
# Objective: minimise sum of x[r][w] for all rounds and words.
#
# Usage:
#   pip install pulp
#   python3 milp_krakken.py
#
# --------------------------------------------------

from pulp import *

# ==================================================
# CONFIGURATION
# ==================================================

ROUNDS       = 10    # number of permutation rounds to model
NUM_WORDS    = 32    # Krakken state size in 64-bit words
BRANCH_NUM   = 9     # MDS branch number (8x8 MDS over GF(2^8) → BN = 9)
NUM_MDS_ROWS = 4     # 4 MDS rows (one per y-coordinate)
NUM_COLS     = 8     # 8 theta/ARX columns

# ==================================================
# PRECOMPUTED STRUCTURAL MAPS
# ==================================================

# --- Pi permutation ---
# pi maps state word at position i  →  new position pi_out[i]
# i = x*4+y,  new_x = (x + 3*y) % 8,  pi_out[i] = new_x*4 + y
def build_pi():
    pi_out = [0] * NUM_WORDS
    for i in range(NUM_WORDS):
        x = i // 4
        y = i % 4
        new_x = (x + 3 * y) % 8
        pi_out[i] = new_x * 4 + y
    return pi_out

PI_OUT = build_pi()   # PI_OUT[i] = destination of word i after pi

# Inverse pi: PI_IN[j] = source word that ends up at position j after pi
PI_IN = [0] * NUM_WORDS
for i, dst in enumerate(PI_OUT):
    PI_IN[dst] = i

# --- Ink-cloud shuffle ---
# ink_cloud: out[(i*7)%32] = in[i]
# IC_OUT[i] = destination of word i after ink_cloud
IC_OUT = [(i * 7) % 32 for i in range(NUM_WORDS)]

# Inverse: IC_IN[j] = source of output position j after ink_cloud
IC_IN = [0] * NUM_WORDS
for i, dst in enumerate(IC_OUT):
    IC_IN[dst] = i

# --- MDS rows ---
# MDS row y contains words {c*4+y  for c in 0..7}
# These are the words entering the MDS.  After MDS → rho (no position change) → pi,
# these words land at positions {PI_OUT[c*4+y]  for c in 0..7} entering chi of r+1.
MDS_ROWS_IN  = [[c * 4 + y for c in range(8)] for y in range(NUM_MDS_ROWS)]
MDS_ROWS_OUT = [[PI_OUT[c * 4 + y] for c in range(8)] for y in range(NUM_MDS_ROWS)]

# --- Theta columns ---
# Column c = {4c, 4c+1, 4c+2, 4c+3}
THETA_COLS = [[4 * c + y for y in range(4)] for c in range(NUM_COLS)]

# Which MDS row does each word belong to? (= y coordinate = w % 4)
def mds_row_of(w):
    return w % 4

# ==================================================
# MILP MODEL
# ==================================================

model = LpProblem("Krakken_Active_Sboxes", LpMinimize)

def binvar(name):
    return LpVariable(name, lowBound=0, upBound=1, cat=LpBinary)

# ---- S-box input activity: x[r][w] ----
# x[r][w] = 1 if word w is active entering chi (S-box layer) of round r.
x = {}
for r in range(ROUNDS + 1):
    for w in range(NUM_WORDS):
        x[r, w] = binvar(f"x_{r}_{w}")

# ---- S-box output activity: a[r][w] ----
# chi is a bijection → a[r][w] = x[r][w].  We enforce a >= x.
a = {}
for r in range(ROUNDS):
    for w in range(NUM_WORDS):
        a[r, w] = binvar(f"a_{r}_{w}")
        model += a[r, w] >= x[r, w], f"sbox_out_{r}_{w}"

# ---- Butterfly stage activity: bf[r][s][w] ----
# bf[r][s][w] = activity of word w after butterfly stage s (s=0..4) of round r.
# bf[r][-1] ≡ a[r]  (input to butterfly = chi output).
# At each stage s: pairs (i, i^(1<<s)) where (i & (1<<s)) == 0 are XOR-mixed.
# Output of each word in the pair is active if either input is active.
bf = {}
for r in range(ROUNDS):
    for s in range(5):
        for w in range(NUM_WORDS):
            bf[r, s, w] = binvar(f"bf_{r}_{s}_{w}")

for r in range(ROUNDS):
    for s in range(5):
        dist = 1 << s
        prev = {w: (a[r, w] if s == 0 else bf[r, s - 1, w]) for w in range(NUM_WORDS)}
        for w in range(NUM_WORDS):
            partner = w ^ dist
            # Output w active if input w or input partner active
            model += bf[r, s, w] >= prev[w],       f"bf_self_{r}_{s}_{w}"
            model += bf[r, s, w] >= prev[partner],  f"bf_partner_{r}_{s}_{w}"
            # Output w inactive only if both inputs inactive (upper bound for LP tightness)
            model += bf[r, s, w] <= prev[w] + prev[partner], f"bf_upper_{r}_{s}_{w}"

# After butterfly stage 4 we also go through ARX (mixes within 4-word columns).
# ARX is linear modular arithmetic: conservative model: ARX output word w is active
# if any word in the same 4-word column (same c = w//4) is active after butterfly.
arx = {}
for r in range(ROUNDS):
    for w in range(NUM_WORDS):
        arx[r, w] = binvar(f"arx_{r}_{w}")
    # Column activity indicator after butterfly
    for c in range(NUM_COLS):
        col_words = THETA_COLS[c]  # same grouping as ARX columns
        for w in col_words:
            # word w after ARX active if any col word was active after butterfly
            for w2 in col_words:
                model += arx[r, w] >= bf[r, 4, w2], f"arx_spread_{r}_{c}_{w}_{w2}"
            # ARX cannot activate words that weren't active (upper bound)
        # Collectively: if any col word active after bf, all may activate;
        # if none active, none activate → sum(arx col) <= sum(bf col) * |col|
        model += lpSum(arx[r, w] for w in col_words) <= \
                 4 * lpSum(bf[r, 4, w] for w in col_words), \
                 f"arx_bound_{r}_{c}"

# ---- Ink-cloud activity: t[r][w] ----
# ink_cloud: t[r][(i*7)%32] = arx[r][i]
t = {}
for r in range(ROUNDS):
    for w in range(NUM_WORDS):
        t[r, w] = binvar(f"t_{r}_{w}")
    for i in range(NUM_WORDS):
        dst = IC_OUT[i]
        model += t[r, dst] == arx[r, i], f"inkcloud_{r}_{i}"

# ---- Theta column indicators: uc[r][c] ----
# uc[r][c] = 1 if any word in column c is active after ink_cloud.
uc = {}
for r in range(ROUNDS):
    for c in range(NUM_COLS):
        uc[r, c] = binvar(f"uc_{r}_{c}")
        col_words = THETA_COLS[c]
        # uc = 1 if any word active
        for w in col_words:
            model += uc[r, c] >= t[r, w], f"uc_lower_{r}_{c}_{w}"
        model += uc[r, c] <= lpSum(t[r, w] for w in col_words), f"uc_upper_{r}_{c}"

# ---- Post-theta activity: v[r][w] ----
# theta_scalar: d[c] = rotr(parity[c-1], 1) XOR parity[c+1]
#               state[4c+y] ^= d[c]
# So column c can be activated by columns (c-1)%8 and (c+1)%8.
# v[r][w] >= t[r][w]                           (own activity preserved)
# v[r][4c+y] >= uc[r][(c-1)%8]                (spread from left neighbour)
# v[r][4c+y] >= uc[r][(c+1)%8]                (spread from right neighbour)
v = {}
for r in range(ROUNDS):
    for w in range(NUM_WORDS):
        v[r, w] = binvar(f"v_{r}_{w}")
    for c in range(NUM_COLS):
        for y in range(4):
            w = 4 * c + y
            model += v[r, w] >= t[r, w],               f"theta_self_{r}_{c}_{y}"
            model += v[r, w] >= uc[r, (c - 1) % 8],    f"theta_left_{r}_{c}_{y}"
            model += v[r, w] >= uc[r, (c + 1) % 8],    f"theta_right_{r}_{c}_{y}"

# ---- MDS row activity indicators: mds_act[r][row] ----
# The MDS (tentacle_mds) operates on each of the 4 rows of 8 words.
# Branch-number constraint: sum(active inputs) + sum(active outputs) >= 9 * indicator
# MDS input  for row y: v[r][c*4+y]         for c in 0..7  (= MDS_ROWS_IN[y])
# MDS output for row y: x[r+1][PI_OUT[c*4+y]] for c in 0..7  (= MDS_ROWS_OUT[y])
mds_act = {}
for r in range(ROUNDS):
    for row in range(NUM_MDS_ROWS):
        mds_act[r, row] = binvar(f"mds_{r}_{row}")
        in_words  = MDS_ROWS_IN[row]
        out_words = MDS_ROWS_OUT[row]
        # Branch-number constraint
        model += (lpSum(v[r, w] for w in in_words) +
                  lpSum(x[r + 1, w] for w in out_words) >=
                  BRANCH_NUM * mds_act[r, row]), f"mds_bn_{r}_{row}"
        # Row active if any input or output word is active
        for w in in_words:
            model += mds_act[r, row] >= v[r, w], f"mds_in_act_{r}_{row}_{w}"
        for w in out_words:
            model += mds_act[r, row] >= x[r + 1, w], f"mds_out_act_{r}_{row}_{w}"

# ==================================================
# OBJECTIVE
# ==================================================
# Minimise total active S-box words over all rounds (= S-box input activity).
# Each active word contributes up to 8 active S-boxes (byte-level); this model
# finds the minimum number of active words, giving a lower bound on active S-boxes.
# For a byte-level count multiply the result by 1 (word active ↔ ≥1 byte active).
model += lpSum(x[r, w] for r in range(ROUNDS) for w in range(NUM_WORDS))

# ==================================================
# NON-ZERO INPUT CONSTRAINT
# ==================================================
model += lpSum(x[0, w] for w in range(NUM_WORDS)) >= 1

# ==================================================
# SOLVE
# ==================================================

print("=" * 60)
print("MILP ACTIVE S-BOX SEARCH  —  Krakken-2048")
print("=" * 60)
print(f"Rounds      : {ROUNDS}")
print(f"Words/round : {NUM_WORDS}  ({NUM_WORDS * 8} S-boxes at byte level)")
print(f"MDS BN      : {BRANCH_NUM}")
print()

solver = PULP_CBC_CMD(msg=True)
model.solve(solver)

print()
print("Status:", LpStatus[model.status])

if model.status != 1:
    print("No feasible solution found.")
    exit(1)

# ==================================================
# RESULTS
# ==================================================

total_active_words = 0
print("\nPer-round active words (word = 64-bit / 8 S-boxes):\n")

for r in range(ROUNDS):
    active_w = sum(1 for w in range(NUM_WORDS) if value(x[r, w]) > 0.5)
    total_active_words += active_w
    # Minimum active S-boxes from this round = active_w (≥1 byte per active word).
    # If MDS guarantees full word activation, count can be 8×active_w in the
    # worst case, but the MILP itself gives the exact lower bound at word level.
    print(f"  Round {r + 1:2d} : {active_w:2d} active words  "
          f"(>= {active_w} S-boxes, up to {active_w * 8} S-boxes)")

print()
print("=" * 60)
print(f"MINIMUM ACTIVE WORDS  = {total_active_words}")
print(f"  → lower bound on active S-boxes: {total_active_words}")
print(f"  → upper bound (all bytes active): {total_active_words * 8}")
print("=" * 60)

# ==================================================
# SHOW TRAIL (word level)
# ==================================================

print("\nMinimal activity trail (word level, 1=active):\n")
header = "Round |  " + "  ".join(f"W{w:02d}" for w in range(NUM_WORDS))
print(header)
print("-" * len(header))
for r in range(ROUNDS + 1):
    row_vals = [int(value(x[r, w]) + 0.5) for w in range(NUM_WORDS)]
    row_str = "  ".join(str(v) for v in row_vals)
    print(f"  {r:3d}  |  {row_str}")

print()
print("=" * 60)
