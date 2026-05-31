import pulp

def solve_milp():
    # 2 Rounds of Krakken S-box search
    prob = pulp.LpProblem("Krakken_Active_Sboxes", pulp.LpMinimize)

    # State variables: for each round r (0, 1), word w (0..31), byte b (0..7)
    # x[r][w][b] is the activity before Chi (non-linear S-box layer)
    # y[r][w][b] is the activity after Chi
    x = {}
    y = {}
    for r in range(2):
        x[r] = {}
        y[r] = {}
        for w in range(32):
            x[r][w] = [pulp.LpVariable(f"x_{r}_{w}_{b}", cat="Binary") for b in range(8)]
            y[r][w] = [pulp.LpVariable(f"y_{r}_{w}_{b}", cat="Binary") for b in range(8)]

    # Chi layer constraints:
    # Pairs (2p, 2p+1)
    # a' = S(a ^ rotl(b, 32))
    # b' = S(b ^ rotl(a', 32))
    # S-box is active if input is active.
    # So y[r][2p][b] = x[r][2p][b] OR x[r][2p+1][(b+4)%8]
    # y[r][2p+1][b] = x[r][2p+1][b] OR y[r][2p][(b+4)%8]
    active_sboxes = []
    for r in range(2):
        for p in range(16):
            w_a = 2 * p
            w_b = 2 * p + 1
            for b in range(8):
                # S-box 1 (for a')
                s1 = pulp.LpVariable(f"s1_{r}_{p}_{b}", cat="Binary")
                prob += s1 >= x[r][w_a][b]
                prob += s1 >= x[r][w_b][(b+4)%8]
                prob += s1 <= x[r][w_a][b] + x[r][w_b][(b+4)%8]
                prob += y[r][w_a][b] == s1
                active_sboxes.append(s1)

                # S-box 2 (for b')
                s2 = pulp.LpVariable(f"s2_{r}_{p}_{b}", cat="Binary")
                prob += s2 >= x[r][w_b][b]
                prob += s2 >= y[r][w_a][(b+4)%8]
                prob += s2 <= x[r][w_b][b] + y[r][w_a][(b+4)%8]
                prob += y[r][w_b][b] == s2
                active_sboxes.append(s2)

    # Butterfly layer + MDS layer between round 0 and round 1:
    # We model MDS. MDS operates per row y (0..3) and per byte b (0..8) on 8 words.
    # The words in row y are 4*c + y for c in 0..7.
    # Input to MDS is the state after Pi/Butterfly/etc.
    # For a simple, conservative, and fast model, we can map the output of round 0 Chi (y[0])
    # to the input of round 1 Chi (x[1]) through MDS and permutations.
    # Let's track how bytes map:
    # y[0][w][b] -> permutation -> input to MDS -> MDS -> output of MDS -> input to round 1 Chi.
    # MDS mixes bytes at position b across the 8 columns for each row y.
    # Let mds_in[y][b][c] be the active bytes entering MDS.
    # Let mds_out[y][b][c] be the active bytes leaving MDS.
    # Branch number constraint:
    # sum_c (mds_in + mds_out) >= 9 * active_row_byte
    # mds_in <= active_row_byte, mds_out <= active_row_byte
    mds_in = {}
    mds_out = {}
    active_rb = {}
    for row in range(4):
        mds_in[row] = {}
        mds_out[row] = {}
        active_rb[row] = {}
        for b in range(8):
            mds_in[row][b] = [pulp.LpVariable(f"mi_{row}_{b}_{c}", cat="Binary") for c in range(8)]
            mds_out[row][b] = [pulp.LpVariable(f"mo_{row}_{b}_{c}", cat="Binary") for c in range(8)]
            arb = pulp.LpVariable(f"arb_{row}_{b}", cat="Binary")
            active_rb[row][b] = arb

            # MDS branch number constraints
            prob += pulp.lpSum(mds_in[row][b]) + pulp.lpSum(mds_out[row][b]) >= 9 * arb
            for c in range(8):
                prob += mds_in[row][b][c] <= arb
                prob += mds_out[row][b][c] <= arb

    # Now link y[0][w][b] to mds_in[row][b][c] and mds_out[row][b][c] to x[1][w][b]
    # For simplicity of the model, since Butterfly and Shuffles permute words and bytes,
    # we can map y[0][w][b] directly to mds_in[row][b][c] via the permutation layers.
    # Specifically:
    # 1. After Chi: y[0][w][b]
    # 2. Butterfly Diffusion: XORs and mixes words. For a differential model,
    #    an active byte at y[0][w][b] can propagate to its butterfly paired word.
    #    To make the model simple and fast, we map each output of MDS to the next round,
    #    and link the activities.
    #    Let's assume the linear layers (Theta, Butterfly, Shuffles) distribute the differences.
    #    A standard conservative model is to connect them directly:
    #    mds_in[row][b][c] corresponds to some byte of word w before MDS.
    #    Let's link:
    #    mds_in[w % 4][b][w // 4] == y[0][w][b]
    #    and
    #    x[1][w][b] == mds_out[w % 4][b][w // 4]
    #    (This represents the row-column mapping: w = 4*c + y, so row = w % 4, c = w // 4).
    for w in range(32):
        row = w % 4
        c = w // 4
        for b in range(8):
            prob += mds_in[row][b][c] == y[0][w][b]
            prob += x[1][w][b] == mds_out[row][b][c]

    # Objective: Minimize active S-boxes in round 0 and round 1
    prob += pulp.lpSum(active_sboxes)

    # Constraint: At least one S-box is active (non-zero difference)
    prob += pulp.lpSum(active_sboxes) >= 1

    # Solve
    prob.solve(pulp.PULP_CBC_CMD(msg=False))

    print("=== MILP S-box Search Results ===")
    print(f"Status: {pulp.LpStatus[prob.status]}")
    print(f"Minimum active S-boxes over 2 rounds: {int(pulp.value(prob.objective))}")

if __name__ == "__main__":
    solve_milp()
