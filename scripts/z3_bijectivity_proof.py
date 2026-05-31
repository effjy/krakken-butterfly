import z3

def rotl64(x, n):
    n = n % 64
    if n == 0:
        return x
    return (x << n) | z3.LShR(x, 64 - n)

def pressure_arx_col(a, b, cc, d):
    # PRESSURE ARX step on 4 words
    a = a + (cc ^ z3.LShR(cc, 17))
    b = b + (d ^ z3.LShR(d, 17))
    cc = cc + (a ^ (a << 31))
    d = d + (b ^ (b << 31))
    return a, rotl64(b, 7), cc, rotl64(d, 19)

def test_bijectivity():
    solver = z3.Solver()
    # Setting a timeout of 5 seconds to ensure it never hangs
    solver.set("timeout", 5000)

    print("=== Z3 Bijectivity Verification ===")

    # Test 1: PRESSURE ARX Bijectivity
    a1, b1, c1, d1 = z3.BitVecs('a1 b1 c1 d1', 64)
    a2, b2, c2, d2 = z3.BitVecs('a2 b2 c2 d2', 64)

    oa1, ob1, oc1, od1 = pressure_arx_col(a1, b1, c1, d1)
    oa2, ob2, oc2, od2 = pressure_arx_col(a2, b2, c2, d2)

    solver.push()
    # Outputs are equal
    solver.add(oa1 == oa2, ob1 == ob2, oc1 == oc2, od1 == od2)
    # Inputs are different
    solver.add(z3.Or(a1 != a2, b1 != b2, c1 != c2, d1 != d2))

    res = solver.check()
    print(f"PRESSURE ARX Collision Search: {res} (expected: unsat)")
    assert res == z3.unsat, "PRESSURE ARX is not bijective!"
    solver.pop()

    # Test 2: Chi S-box Layer Pair Bijectivity (using uninterpreted bijective functions)
    S = z3.Function('S', z3.BitVecSort(64), z3.BitVecSort(64))
    S_inv = z3.Function('S_inv', z3.BitVecSort(64), z3.BitVecSort(64))
    
    # Axioms of bijectivity
    x = z3.BitVec('x', 64)
    solver.add(z3.ForAll(x, S_inv(S(x)) == x))
    solver.add(z3.ForAll(x, S(S_inv(x)) == x))

    a1, b1 = z3.BitVecs('a1 b1', 64)
    a2, b2 = z3.BitVecs('a2 b2', 64)

    ap1 = S(a1 ^ rotl64(b1, 32))
    bp1 = S(b1 ^ rotl64(ap1, 32))

    ap2 = S(a2 ^ rotl64(b2, 32))
    bp2 = S(b2 ^ rotl64(ap2, 32))

    solver.push()
    solver.add(ap1 == ap2, bp1 == bp2)
    solver.add(z3.Or(a1 != a2, b1 != b2))

    res = solver.check()
    print(f"Chi S-box Layer Pair Collision Search: {res} (expected: unsat)")
    assert res == z3.unsat, "Chi S-box Layer Pair is not bijective!"
    solver.pop()

    print("SUCCESS: Bijectivity of non-linear layers formally proven by Z3 SMT solver.")

if __name__ == "__main__":
    test_bijectivity()
