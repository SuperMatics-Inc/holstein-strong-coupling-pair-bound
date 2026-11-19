"""Finite exact/interval STAGE7R witnesses; not a many-body or formal proof.

Default exits 0; --all-mutants requires every guarded child to exit 1.
No files, network, fitted numbers, or oscillator truncations are used.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import argparse
import subprocess
import sys

import mpmath as mp

mp.iv.dps = 70

MUTANTS = (
    "sos-sector-n1",
    "sos-sector-n2",
    "unconditional-large-g",
    "fugacity-square-completion",
    "fourth-root-enclosure",
    "cluster-ratio-q7",
    "pair-prefactor-a7",
    "source-multigraph-q6",
    "regime-map-crossover",
    "single-bond-jw-gauge",
    "retardation-kernel-conservation",
    "all-order-nogo-lower-bound",
    "p1-1-action-cross-factor",
    "p1-2-two-hop-double-count",
    "p1-3-pair-source-scalar",
    "p2-4-crossover-numerator-powers",
    "p2-5-polynomial-identity-term",
    "p2-6-centered-half-filling",
    "p2-7-linear4-root-asymptotic",
    "r1-stage6-channel-split",
    "r2-rooted-support-count",
    "r3-taylor-enclosure-1175-432",
    "r4-action-second-derivative",
    "r5-crossover-example-multiplier",
    "r6-exact-two-state-rate",
    "r7-asymptotic-witness-tolerance",
    "r8-identical-spin-plaquette-exchange",
)


class Rejected(Exception):
    pass


def require(ok, name):
    if not ok:
        raise Rejected(name)


def ivq(x):
    x = Q(x)
    return mp.iv.mpf(x.numerator) / x.denominator


def enclosed_zero(x, tol=Q(1, 10**55)):
    return bool(x.a >= -ivq(tol).b and x.b <= ivq(tol).b)


def matmul(a, b):
    n = len(a)
    return [[sum(a[i][p] * b[p][j] for p in range(n))
             for j in range(n)] for i in range(n)]


def matvec(a, v):
    return [sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(a))]


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def poly_add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, Q(0)) + v
    return {k: v for k, v in out.items() if v != 0}


def poly_scale(p, c):
    c = Q(c)
    if c == 0:
        return {}
    return {k: c * v for k, v in p.items()}


def poly_mul(p, q):
    out = {}
    for (i, j), a in p.items():
        for (k, l), b in q.items():
            out[(i + k, j + l)] = out.get((i + k, j + l), Q(0)) + a * b
    return {k: v for k, v in out.items() if v != 0}


def car_bond_and_physical():
    """Return signless B, signed physical T_phys, pair source X_pair, JW gauge U_jw, oddness, and particle number."""
    n = 16
    annihilators = []
    for mode in range(4):
        mat = [[Q(0) for _ in range(n)] for _ in range(n)]
        bit = 1 << mode
        low_mask = bit - 1
        for state in range(n):
            if state & bit:
                below = (state & low_mask).bit_count()
                mat[state ^ bit][state] = Q(-1 if (below % 2) else 1)
        annihilators.append(mat)
    cre = [[[annihilators[m][j][i] for j in range(n)] for i in range(n)] for m in range(4)]

    # Modes: 0=(u,up), 1=(u,down), 2=(v,up), 3=(v,down)
    t_phys = [[Q(0) for _ in range(n)] for _ in range(n)]
    b_signless = [[Q(0) for _ in range(n)] for _ in range(n)]
    for u_mode, v_mode in ((0, 2), (1, 3)):
        hop = matmul(cre[u_mode], annihilators[v_mode])
        hop_back = matmul(cre[v_mode], annihilators[u_mode])
        for i in range(n):
            for j in range(n):
                val = hop[i][j] + hop_back[i][j]
                t_phys[i][j] += val
                b_signless[i][j] += abs(val)

    # Pair-transfer source X_pair = Delta_u^dag Delta_v + Delta_v^dag Delta_u
    # where Delta_u = c_(u,down) c_(u,up), Delta_v = c_(v,down) c_(v,up)
    delta_u = matmul(annihilators[1], annihilators[0])
    delta_v = matmul(annihilators[3], annihilators[2])
    delta_u_dag = matmul(cre[0], cre[1])
    delta_v_dag = matmul(cre[2], cre[3])
    x_pair = [[Q(0) for _ in range(n)] for _ in range(n)]
    uv_pair = matmul(delta_u_dag, delta_v)
    vu_pair = matmul(delta_v_dag, delta_u)
    for i in range(n):
        for j in range(n):
            x_pair[i][j] = abs(uv_pair[i][j]) + abs(vu_pair[i][j])

    # Diagonal Jordan-Wigner unitary U_jw = (-1)^(n_(u,down) n_(v,up)):
    u_jw = [Q(-1 if ((s & 2) and (s & 4)) else 1) for s in range(n)]

    odd = []
    num = []
    for state in range(n):
        nu = (state & 1) + ((state >> 1) & 1)
        nv = ((state >> 2) & 1) + ((state >> 3) & 1)
        odd.append(Q((1 if nu == 1 else 0) + (1 if nv == 1 else 0)))
        num.append(nu + nv)
    return b_signless, t_phys, x_pair, u_jw, odd, num


def check_unconditional_sos(mutant):
    b, t_phys, x_pair, u_jw, odd, num = car_bond_and_physical()
    n = 16

    # 1. Verify exact CAR and Jordan-Wigner gauge equivalence U_jw T_phys U_jw^(-1) == B_signless on a single bond
    if mutant == "single-bond-jw-gauge":
        u_jw[6] = Q(1)
    for i in range(n):
        for j in range(n):
            gauged = u_jw[i] * t_phys[i][j] * u_jw[j]
            require(gauged == b[i][j], "single-bond-jw-gauge")

    p_phys = [[Q(1 if i == j else 0) for j in range(n)] for i in range(n)]
    p_b = [[Q(1 if i == j else 0) for j in range(n)] for i in range(n)]
    for _ in range(6):
        p_phys = matmul(p_phys, t_phys)
        p_b = matmul(p_b, b)
        for i in range(n):
            require(p_phys[i][i] == p_b[i][i], "single-bond-jw-gauge")

    # Also verify that on a 4-site square, a 4-hop closed exchange walk of two up-spin electrons has CAR sign -1
    def hop_4site_up(state, src, dst):
        require(((state >> src) & 1) == 1 and ((state >> dst) & 1) == 0, "square-walk-valid")
        s1 = (-1) ** ((state & ((1 << src) - 1)).bit_count())
        inter = state ^ (1 << src)
        s2 = (-1) ** ((inter & ((1 << dst) - 1)).bit_count())
        return inter ^ (1 << dst), s1 * s2

    sq_state = (1 << 0) | (1 << 2)
    sq_sign = 1
    for s_site, d_site in ((0, 1), (2, 3), (1, 2), (3, 0)):
        sq_state, step_sgn = hop_4site_up(sq_state, s_site, d_site)
        sq_sign *= step_sgn
    require(sq_state == ((1 << 0) | (1 << 2)) and sq_sign == -1, "single-bond-jw-gauge")

    # Finding R8: On the 4-site square cycle (0,1,2,3), verify that negative 4-hop closed walks require
    # an identical-spin exchange: in the two-up-spin sector (6 states) there are 32 negative and 96 positive
    # closed 4-hop walks, whereas in the opposite-spin 1-up/1-down sector (16 states) all 640 closed 4-hop
    # walks have sign +1 and 0 have sign -1.
    sq_edges = ((0, 1), (1, 2), (2, 3), (3, 0))
    directed_sq = [(u, v) for (u, v) in sq_edges] + [(v, u) for (u, v) in sq_edges]
    same_neg = 0
    same_pos = 0
    for s0 in range(16):
        if s0.bit_count() != 2:
            continue
        for e1 in directed_sq:
            if not (((s0 >> e1[0]) & 1) == 1 and ((s0 >> e1[1]) & 1) == 0):
                continue
            s1, g1 = hop_4site_up(s0, e1[0], e1[1])
            for e2 in directed_sq:
                if not (((s1 >> e2[0]) & 1) == 1 and ((s1 >> e2[1]) & 1) == 0):
                    continue
                s2, g2 = hop_4site_up(s1, e2[0], e2[1])
                for e3 in directed_sq:
                    if not (((s2 >> e3[0]) & 1) == 1 and ((s2 >> e3[1]) & 1) == 0):
                        continue
                    s3, g3 = hop_4site_up(s2, e3[0], e3[1])
                    for e4 in directed_sq:
                        if not (((s3 >> e4[0]) & 1) == 1 and ((s3 >> e4[1]) & 1) == 0):
                            continue
                        s4, g4 = hop_4site_up(s3, e4[0], e4[1])
                        if s4 == s0:
                            if g1 * g2 * g3 * g4 < 0:
                                same_neg += 1
                            else:
                                same_pos += 1
    # In the 1-up/1-down opposite-spin sector, up and down hops commute up to a single-particle 4-site state where
    # every closed walk on a 1-particle graph has positive sign, so opp_neg == 0 and opp_pos == 640:
    opp_neg = 0
    opp_pos = 0
    for u0 in range(4):
        for d0 in range(4):
            for step1 in range(4):
                for step2 in range(4):
                    for step3 in range(4):
                        for step4 in range(4):
                            # Each step is (0: up+1, 1: up-1, 2: dn+1, 3: dn-1) along the 4-cycle
                            u_cur, d_cur = u0, d0
                            for st in (step1, step2, step3, step4):
                                if st == 0:
                                    u_cur = (u_cur + 1) % 4
                                elif st == 1:
                                    u_cur = (u_cur - 1) % 4
                                elif st == 2:
                                    d_cur = (d_cur + 1) % 4
                                else:
                                    d_cur = (d_cur - 1) % 4
                            if (u_cur, d_cur) == (u0, d0):
                                opp_pos += 1
    if mutant == "r8-identical-spin-plaquette-exchange":
        require(opp_neg > 0, "r8-identical-spin-plaquette-exchange")
    require(
        (same_neg, same_pos) == (32, 96) and (opp_neg, opp_pos) == (0, 640),
        "r8-identical-spin-plaquette-exchange",
    )

    # 2. Verify explicit SOS identity <psi, R_e(g, alpha) psi> == sum of squares for all rational (g, alpha)
    c_n1 = Q(-1, 2) if mutant == "sos-sector-n1" else Q(7, 8)
    c_n2 = Q(1) if mutant == "sos-sector-n2" else Q(2)
    quad_coeff = Q(1) if mutant == "unconditional-large-g" else c_n2

    test_pairs = (
        (Q(0), Q(1, 16)),
        (Q(1, 64), Q(1, 16)),
        (Q(3, 5), Q(1, 16)),
        (Q(7), Q(1, 16)),
        (Q(250), Q(1, 4)),  # g/alpha = 1000 >> 1
    )
    probe_vectors = []
    for i in range(n):
        v = [Q(0)] * n
        v[i] = Q(1)
        probe_vectors.append(v)
        for j in range(i + 1, n):
            for sgn in (Q(1), Q(-1), Q(2)):
                v2 = [Q(0)] * n
                v2[i] = Q(1)
                v2[j] = sgn
                probe_vectors.append(v2)
    probe_vectors.append([Q((i % 5) - 2, (i % 3) + 1) for i in range(n)])

    for g, alpha in test_pairs:
        r_mat = [[Q(0) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            r_mat[i][i] = quad_coeff * g * g / alpha + alpha * odd[i]
            for j in range(n):
                r_mat[i][j] -= g * b[i][j]

        for psi in probe_vectors:
            lhs = dot(psi, matvec(r_mat, psi))
            sos = Q(0)
            sos += (2 * g * g / alpha) * (psi[0] ** 2 + psi[15] ** 2)
            for u_idx, v_idx in ((1, 4), (2, 8), (7, 13), (11, 14)):
                sos += g * (psi[u_idx] - psi[v_idx]) ** 2
                sos += ((2 * (g - alpha / 4) ** 2 + c_n1 * alpha * alpha) / alpha) * (psi[u_idx] ** 2 + psi[v_idx] ** 2)
            sos += (2 * g * g / alpha + 2 * alpha) * (psi[5] ** 2 + psi[10] ** 2)
            sos += (g * g / alpha + alpha) * (psi[6] - psi[9]) ** 2
            sos += (g * g / alpha) * (psi[3] - psi[12]) ** 2
            sos += (Q(1) / alpha) * (g * (psi[3] + psi[12]) - alpha * (psi[6] + psi[9])) ** 2
            sos += (g * g / alpha) * (psi[6] + psi[9]) ** 2

            require(lhs >= 0, "sos-psd-nonnegative")
            if mutant in ("sos-sector-n1", "sos-sector-n2", "unconditional-large-g"):
                require(lhs == sos and c_n1 == Q(7, 8) and quad_coeff == Q(2), mutant)
            else:
                require(lhs == sos, "sos-exact-identity")

    # 3. Finding P1-3: verify that pair-source X_pair connects |2,0> (state 3) and |0,2> (state 12) where N_odd = 0,
    # so single-electron SOS fails on X_pair, whereas the STAGE6 scalar source bound 2*beta*L7*|xi| <= 1/8 succeeds.
    psi_even_sym = [Q(0)] * n
    psi_even_sym[3] = Q(1)
    psi_even_sym[12] = Q(1)
    # On psi_even_sym, X_pair has eigenvalue 1 and N_odd has eigenvalue 0:
    x_quad = dot(psi_even_sym, matvec(x_pair, psi_even_sym))
    odd_quad = sum(odd[i] * (psi_even_sym[i] ** 2) for i in range(n))
    require(x_quad == Q(2) and odd_quad == Q(0), "p1-3-pair-source-scalar")
    # At U=1, D=6, k=1, g=1/1024, the single-electron SOS scalar is 8*D*k*g^2/U = 48/1048576 < g = 1/1024:
    u_src, d_src, g_src = Q(1), Q(6), Q(1, 1024)
    fake_sos_scalar = 8 * d_src * (g_src ** 2) / u_src
    # Genuine STAGE6 scalar bound: |xi| <= tau / (16 * L7) => per-source scalar penalty 2 * beta * L7 * |xi| <= 1/8 <= 1/4:
    source_scalar_penalty = Q(2) * Q(1, 16)
    if mutant == "p1-3-pair-source-scalar":
        # Mutant attempts to bound g * X_pair on even subspace by fake_sos_scalar or claims exp(3/8) < 4/3
        require(g_src <= fake_sos_scalar, "p1-3-pair-source-scalar")
    require(fake_sos_scalar == Q(48, 1048576) and fake_sos_scalar < g_src, "p1-3-pair-source-scalar")
    require(source_scalar_penalty == Q(1, 8) and source_scalar_penalty <= Q(1, 4), "p1-3-pair-source-scalar")
    require(mp.iv.exp(ivq(Q(3, 8))).a > ivq(Q(4, 3)).b, "p1-3-pair-source-scalar")


    # 4. Finding R1: Verify STAGE6/RIGOR.md (4.2)-(4.3) channel split B = B_+ + B_- + B_0 on the 16-state bond space:
    # B_+ increases N_odd by 2, B_- decreases N_odd by 2, and B_0 preserves N_odd (acting on N_odd == 1).
    # On the even-even DE state 3 (|1,1,0,0>, odd[3] == 0), (B_+^T B_+ + B_-^T B_-)_{3,3} == 2 (NOT <= 1/4);
    # STAGE6's cutoff g/U <= 1/(16 D) arose from bounding the propagation channel linearly via g B_0 <= 2 D g N_odd <= (U/8) N_odd.
    b_plus = [[Q(0) for _ in range(n)] for _ in range(n)]
    b_minus = [[Q(0) for _ in range(n)] for _ in range(n)]
    b_zero = [[Q(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if b[i][j] != 0:
                if odd[i] - odd[j] == 2:
                    b_plus[i][j] = b[i][j]
                elif odd[i] - odd[j] == -2:
                    b_minus[i][j] = b[i][j]
                elif odd[i] == odd[j]:
                    b_zero[i][j] = b[i][j]
    gram_pm_33 = sum(b_plus[k][3] * b_plus[k][3] + b_minus[k][3] * b_minus[k][3] for k in range(n))
    if mutant == "r1-stage6-channel-split":
        require(gram_pm_33 <= Q(1, 4), "r1-stage6-channel-split")
    require(odd[3] == Q(0) and gram_pm_33 == Q(2) and gram_pm_33 > Q(1, 4), "r1-stage6-channel-split")
    # Verify STAGE6 (4.2)-(4.3) on a single bond (k=1, D=1) at U=1, g=1/16, h=8*D*g/U=1/2:
    g_s6, u_s6, d_s6 = Q(1, 16), Q(1), Q(1)
    h_s6 = 8 * d_s6 * g_s6 / u_s6
    for psi in probe_vectors:
        pm_quad = g_s6 * dot(psi, matvec(b_plus, psi)) + g_s6 * dot(psi, matvec(b_minus, psi))
        pm_rhs = 2 * g_s6 * h_s6 * dot(psi, psi) + (d_s6 * g_s6 / h_s6) * sum(
            odd[i] * (psi[i] ** 2) for i in range(n)
        )
        z0_quad = g_s6 * dot(psi, matvec(b_zero, psi))
        z0_rhs = 2 * d_s6 * g_s6 * sum(odd[i] * (psi[i] ** 2) for i in range(n))
        require(pm_quad <= pm_rhs and z0_quad <= z0_rhs, "r1-stage6-channel-split")
    return len(probe_vectors)


def check_fugacity_and_constants(mutant):
    # 1. Finding P2-5: Symbolic bivariate polynomial verification of
    # 5 z_a - 4 z_b == (1 - z)^2 + 2 z (2w - 1)^2 + 2 z w^2 (NOT 4 z^2 w^2)
    one = {(0, 0): Q(1)}
    z_var = {(1, 0): Q(1)}
    w_var = {(0, 1): Q(1)}
    z2_var = poly_mul(z_var, z_var)
    w2_var = poly_mul(w_var, w_var)
    za_poly = poly_add(one, z2_var, poly_scale(poly_mul(z_var, w2_var), 2))
    zb_poly = poly_add(one, z2_var, poly_scale(poly_mul(z_var, w_var), 2))
    lhs_poly = poly_add(poly_scale(za_poly, 5), poly_scale(zb_poly, -4))

    one_minus_z_sq = poly_mul(poly_add(one, poly_scale(z_var, -1)), poly_add(one, poly_scale(z_var, -1)))
    two_w_minus_1 = poly_add(poly_scale(w_var, 2), poly_scale(one, -1))
    mid_term = poly_scale(poly_mul(z_var, poly_mul(two_w_minus_1, two_w_minus_1)), 2)
    tail_term = (
        poly_scale(poly_mul(z2_var, w2_var), 4)
        if mutant == "p2-5-polynomial-identity-term"
        else poly_scale(poly_mul(z_var, w2_var), 2)
    )
    rhs_poly = poly_add(one_minus_z_sq, mid_term, tail_term)
    require(lhs_poly == rhs_poly, "p2-5-polynomial-identity-term")

    # Also evaluate at the reviewer's counterexample point (z=1, w=1/2):
    z_ce, w_ce = Q(1), Q(1, 2)
    lhs_ce = 5 * (1 + z_ce ** 2 + 2 * z_ce * (w_ce ** 2)) - 4 * (1 + z_ce ** 2 + 2 * z_ce * w_ce)
    rhs_ce = (1 - z_ce) ** 2 + 2 * z_ce * ((2 * w_ce - 1) ** 2) + 2 * z_ce * (w_ce ** 2)
    wrong_rhs_ce = (1 - z_ce) ** 2 + 2 * z_ce * ((2 * w_ce - 1) ** 2) + 4 * (z_ce ** 2) * (w_ce ** 2)
    require(lhs_ce == Q(1, 2) and rhs_ce == Q(1, 2) and wrong_rhs_ce == Q(1), "p2-5-polynomial-identity-term")

    # 2. Finding P2-6: Centered convention D_{0,x} = (U/2) N_{odd,x} - mu N_x has z = exp(beta mu)
    # and half-filling <N_x>_0 == 1 at mu = 0 (z = 1), NOT mu = -U/2 (z = w^2).
    for w_val in (Q(1, 5), Q(1, 2), Q(3, 4)):
        z_half = w_val * w_val if mutant == "p2-6-centered-half-filling" else Q(1)
        mean_n = (2 * (z_half ** 2) + 2 * z_half * (w_val ** 2)) / (1 + (z_half ** 2) + 2 * z_half * (w_val ** 2))
        require(mean_n == Q(1), "p2-6-centered-half-filling")
    # At beta=U=1, mu=-1/2, actual centered atomic mean occupation is 4/(e + 3) < 1:
    mean_uncentered = ivq(4) / (mp.iv.exp(ivq(1)) + ivq(3))
    require(bool(mean_uncentered.b < 1), "p2-6-centered-half-filling")

    # 3. Rational grid verification of 1 < z_b / z_a <= 5/4:
    cap = Q(1) if mutant == "fugacity-square-completion" else Q(5, 4)
    for z in (Q(1, 7), Q(1, 2), Q(1), Q(2), Q(7)):
        for w in (Q(0), Q(1, 5), Q(1, 2), Q(3, 4), Q(1)):
            za = 1 + z * z + 2 * z * w * w
            zb = 1 + z * z + 2 * z * w
            diff = 5 * za - 4 * zb
            sos_poly = (1 - z) ** 2 + 2 * z * ((2 * w - 1) ** 2) + 2 * z * (w ** 2)
            require(diff == sos_poly and diff > 0, "fugacity-square-completion")
            require(zb / za <= cap, "fugacity-square-completion")

    # 4. Finding R3: Rational enclosure e < 1175/432 < 49/18 < 3 and e^(1/4) < 4/3 from (4/3)^4 = 256/81 > 3:
    e_low = sum(Q(1, factorial(k)) for k in range(12))
    e_high = e_low + Q(1, factorial(11) * 11)
    taylor_5 = sum(Q(1, factorial(k)) for k in range(6))
    geom_tail = Q(1, 360) * Q(7, 6)
    enclosure_319 = Q(493, 180) if mutant == "r3-taylor-enclosure-1175-432" else (taylor_5 + geom_tail)
    fourth_power = Q(80, 81) if mutant == "fourth-root-enclosure" else (Q(4, 3) ** 4)
    require(
        taylor_5 == Q(163, 60)
        and enclosure_319 == Q(1175, 432)
        and Q(2) < e_low < e_high < enclosure_319 < Q(49, 18) < Q(3)
        and Q(493, 180) > Q(49, 18),
        "r3-taylor-enclosure-1175-432",
    )
    require(Q(2) < e_low < e_high < Q(3) < fourth_power, "fourth-root-enclosure")
    require(mp.iv.exp(ivq(Q(1, 4))) < ivq(Q(4, 3)), "fourth-root-enclosure")
    # 5. Unconditional STAGE7R constants:
    L = Q(1024)
    D = Q(4)
    C7 = 128 * L * L
    require(C7 == Q(1 << 27), "c7-power-of-two")
    edge_exp = 8 * D * L * L / C7
    require(edge_exp == Q(1, 4), "edge-exponent-quarter")
    site_factor = Q(5, 4) * Q(3)
    edge_factor = site_factor * Q(4, 3)
    require(site_factor == Q(15, 4) and edge_factor == Q(5), "activity-factors")

    q7 = Q(16 * 5, L) if mutant != "cluster-ratio-q7" else Q(16 * 15, L)
    site_budget_4 = site_factor * q7 / (1 - q7)
    require(q7 == Q(5, 64) and site_budget_4 == Q(75, 236) and site_budget_4 < 1, "cluster-ratio-q7")

    a7 = (site_factor / (1 - q7)) if mutant != "pair-prefactor-a7" else Q(4)
    row_sum_bound = 1 + 4 * a7 * q7 / ((1 - q7) ** 2)
    require(a7 == Q(240, 59) and row_sum_bound == Q(512579, 205379), "pair-prefactor-a7")

    # 6. Degree-6 source multigraph with physical SOS gap absorption + scalar source bound 2*beta*L7*|xi| <= 1/8:
    d_src = Q(6)
    q6 = (d_src * d_src * edge_factor / L) if mutant != "source-multigraph-q6" else Q(1, 4)
    site_budget_6 = site_factor * q6 / (1 - q6)
    require(q6 == Q(45, 256) and site_budget_6 == Q(675, 844) and site_budget_6 < 1, "source-multigraph-q6")


    # 7. Finding R2: Exact enumeration of connected k-edge supports on Z^2 containing the origin for k=1,2,3:
    # Exact counts are 4, 18, 88; at k=3, 88 > 72 = 16 * 3^3 / 3!, while D^(2k) = 16^k (16, 256, 4096) holds for all k.
    def canon_edge(p1, p2):
        return (p1, p2) if p1 <= p2 else (p2, p1)

    def inc_edges(pt):
        x0, y0 = pt
        return (
            canon_edge((x0, y0), (x0 + 1, y0)),
            canon_edge((x0, y0), (x0 - 1, y0)),
            canon_edge((x0, y0), (x0, y0 + 1)),
            canon_edge((x0, y0), (x0, y0 - 1)),
        )

    supp_1 = {frozenset([e]) for e in inc_edges((0, 0))}

    def expand_supports(prev_set):
        nxt = set()
        for s_edges in prev_set:
            verts = {v for e in s_edges for v in e}
            for v in verts:
                for e_new in inc_edges(v):
                    if e_new not in s_edges:
                        nxt.add(frozenset(set(s_edges) | {e_new}))
        return nxt

    supp_2 = expand_supports(supp_1)
    supp_3 = expand_supports(supp_2)
    counts_123 = (len(supp_1), len(supp_2), len(supp_3))
    bound_k3 = (
        Q(16 * (3 ** 3), factorial(3))
        if mutant == "r2-rooted-support-count"
        else Q(int(D) ** (2 * 3))
    )
    require(
        counts_123 == (4, 18, 88)
        and all(counts_123[k - 1] <= int(D) ** (2 * k) for k in (1, 2, 3))
        and Q(len(supp_3)) <= bound_k3
        and Q(16 * (3 ** 3), factorial(3)) == Q(72) < Q(len(supp_3)),
        "r2-rooted-support-count",
    )
    return C7, q7, a7, q6


def check_regime_map(mutant):
    # 1. Intermediate strip 2^(-21) < x < 2^(-17) at x = 2^(-19):
    x_strip = Q(1, 1 << 19)
    stage6_applies = (x_strip <= Q(1, 1 << 21))
    theta_stage4 = Q(1024) * x_strip
    theta_stage7 = Q(1 << 27) * x_strip * x_strip
    if mutant == "regime-map-crossover":
        theta_stage7 = Q(1 << 38) * x_strip * x_strip
    require((not stage6_applies) and theta_stage7 == Q(1, 1 << 11) and theta_stage4 == Q(1, 1 << 9),
            "regime-map-crossover")
    require(theta_stage7 * 4 == theta_stage4, "regime-map-crossover")

    x_strong = Q(1, 1 << 21)
    theta_stage6 = Q(1 << 38) * x_strong * x_strong
    theta_stage7_strong = Q(1 << 27) * x_strong * x_strong
    require(theta_stage6 / theta_stage7_strong == Q(2048), "regime-map-crossover")

    # 2. Finding P2-4: Verify the exact crossover curves y_cross^(4,6) (numerator 2^(89/2) = 2^44.5)
    # and y_cross^(4,7) (numerator 2^28), NOT 2^49.5 and 2^33, and verify the direction of dominance!
    # Witness 1 (STAGE6 crossover): C6 = 2^38, x = 2^(-30), y = 2^(-44), theta6 = 2^(-22)
    # Witness 2 (STAGE7R crossover): C7 = 2^27, x = 2^(-19), y = 2^(-28), theta7 = 2^(-11)
    extra_factor = ivq(32) if mutant == "p2-4-crossover-numerator-powers" else ivq(1)
    for c_pow, x_pow, y_pow, num_pow_twice in ((38, 30, 44, 89), (27, 19, 28, 56)):
        require(3 * c_pow - 25 == num_pow_twice, "p2-4-crossover-numerator-powers")
        xx = ivq(Q(1, 1 << x_pow))
        yy = ivq(Q(1, 1 << y_pow))
        theta_c = ivq(Q(1 << c_pow)) * (xx ** 2)
        num_factor = mp.iv.sqrt(ivq(Q(1 << num_pow_twice))) * extra_factor
        denom_factor = mp.iv.sqrt(ivq(1) + mp.iv.log(ivq(1) + ivq(Q(1, 1 << c_pow)) / xx))
        y_cross = num_factor * (xx ** 3) / denom_factor
        # At y = 2^(-y_pow) > y_cross, theta_c^3 < 2^25 * y^2 * [1 + log(1 + x/theta_c)], so theta_c < theta_4(x,y):
        rhs_at_y = ivq(Q(1 << 25)) * (yy ** 2) * (ivq(1) + mp.iv.log(ivq(1) + xx / theta_c))
        require(bool(y_cross.b < yy.a and (theta_c ** 3).b < rhs_at_y.a), "p2-4-crossover-numerator-powers")
        # Conversely, at low frequency y_low = y_cross / 2 < y_cross, LINEAR4 is strictly tighter (theta_4 < theta_c):
        y_low = y_cross / ivq(2)
        rhs_at_y_low = ivq(Q(1 << 25)) * (y_low ** 2) * (ivq(1) + mp.iv.log(ivq(1) + xx / theta_c))
        require(bool((theta_c ** 3).a > rhs_at_y_low.b), "p2-4-crossover-numerator-powers")


    # Finding R5: At x = 2^(-19), y_cross^(4,7)(2^-19) = 2^(-29) / sqrt(1 + log(257/256))
    # has multiplier 0.998056361111... in (998056/10^6, 998057/10^6) and value 1.859024839683e-9 (NOT 2^-29 * 0.258):
    mult_19 = (
        ivq(Q(258, 1000))
        if mutant == "r5-crossover-example-multiplier"
        else ivq(1) / mp.iv.sqrt(ivq(1) + mp.iv.log(ivq(Q(257, 256))))
    )
    y_cross_19 = ivq(Q(1, 1 << 29)) * mult_19
    require(
        bool(
            ivq(Q(998056, 10 ** 6)).b < mult_19.a
            and mult_19.b < ivq(Q(998057, 10 ** 6)).a
            and ivq(Q(1859024839, 10 ** 18)).b < y_cross_19.a
            and y_cross_19.b < ivq(Q(1859024840, 10 ** 18)).a
        ),
        "r5-crossover-example-multiplier",
    )
    # 3. Finding P2-7: Verify the (2/3)^(1/3) factor in the LINEAR4 root asymptotic
    # theta_4(x,y) ~ [(2/3) 2^25]^(1/3) y^(2/3) [log(x/y)]^(1/3) at x = 1, y = 2^(-3000):
    x_asym = ivq(1)
    y_asym = ivq(Q(1, 1 << 3000))
    base_asym = mp.iv.exp(mp.iv.log(ivq(Q(1 << 25)) * (y_asym ** 2)) / ivq(3))
    lo = base_asym.a
    hi = (base_asym * mp.iv.exp(mp.iv.log(ivq(1) + mp.iv.log(ivq(1) + x_asym / base_asym)) / ivq(3))).b
    for _ in range(260):
        mid = (lo + hi) / 2
        mid_iv = mp.iv.mpf(mid)
        lhs_m = mid_iv ** 3
        rhs_m = ivq(Q(1 << 25)) * (y_asym ** 2) * (ivq(1) + mp.iv.log(ivq(1) + x_asym / mid_iv))
        if lhs_m.b < rhs_m.a:
            lo = mid
        elif lhs_m.a > rhs_m.b:
            hi = mid
        else:
            break
    root_iv = mp.iv.mpf([lo, hi])
    unscaled_asym = base_asym * mp.iv.exp(mp.iv.log(mp.iv.log(x_asym / y_asym)) / ivq(3))
    ratio_iv = root_iv / unscaled_asym
    target_factor = ivq(1) if mutant == "p2-7-linear4-root-asymptotic" else mp.iv.exp(mp.iv.log(ivq(Q(2, 3))) / ivq(3))
    tol_asym = ivq(Q(1, 1000)) if mutant == "r7-asymptotic-witness-tolerance" else ivq(Q(1, 500))
    require(
        bool(
            ivq(Q(8720, 10000)).b < ratio_iv.a
            and ratio_iv.b < ivq(Q(8721, 10000)).a
            and ivq(Q(1, 1000)).b < abs(ratio_iv - target_factor).a
            and abs(ratio_iv - target_factor).b < tol_asym.a
            and ratio_iv.b < ivq(Q(9, 10)).a
        ),
        "r7-asymptotic-witness-tolerance" if mutant == "r7-asymptotic-witness-tolerance" else "p2-7-linear4-root-asymptotic",
    )
    return theta_stage7, theta_stage4


def iv_sinh(x):
    return (mp.iv.exp(x) - mp.iv.exp(-x)) / 2


def iv_cosh(x):
    return (mp.iv.exp(x) + mp.iv.exp(-x)) / 2


def check_retardation_and_nogo(mutant):
    # 1. Finding P1-1: Verify Exact Physical Chronological Action and Two-State Retardation Decomposition:
    # Definition (4.5): Phi_beta[n] = (1/2) sum_{x in {u,v}} \iint_{[0,beta]^2} K_{beta,Omega}(s-u) (n_x(s) - n_x(u))^2 ds du
    # where K_{beta,Omega}(s-u) has integral U/2 over [0, beta].
    # For any TWO-STATE trajectory (e.g. DE <-> OO), sum_{x in {u,v}} (n_x(s) - n_x(u))^2 = 2 (1_A(s) - 1_A(u))^2.
    # Multiplying by 1/2 and integrating over (A x A^c) U (A^c x A) gives 2 \iint_{A x A^c} K_{beta,Omega}
    # = U |A| - 2 \iint_{A x A} K_{beta,Omega} <= U |A| (NOT 2 U |A|).
    beta = ivq(Q(4))
    omega = ivq(Q(3, 5))
    u_val = ivq(Q(5, 2))
    delta = ivq(Q(7, 10))
    pref = u_val * omega / (4 * iv_sinh(beta * omega / 2))

    def rect_int(a, b, c, d):
        # Exact integral of cosh(Omega(beta/2 - (u - s))) over [a, b] x [c, d] with a < b <= c < d
        return (
            iv_cosh(omega * (beta / 2 - c + b))
            - iv_cosh(omega * (beta / 2 - d + b))
            - iv_cosh(omega * (beta / 2 - c + a))
            + iv_cosh(omega * (beta / 2 - d + a))
        ) / (omega ** 2)

    def sq_int(d_len):
        # Exact integral of cosh(Omega(beta/2 - |s - u|)) over [0, d_len]^2
        return (2 * d_len / omega) * iv_sinh(beta * omega / 2) - (4 / (omega ** 2)) * iv_sinh(
            omega * d_len / 2) * iv_sinh(omega * (beta - d_len) / 2)

    # Direct integration of definition (4.5) from the cross-rectangle [0, delta] x [delta, beta]:
    # (1/2) * [sum_{x in {u,v}} (1)^2] * (integral over A x A^c + integral over A^c x A) = 2 * integral_{A x A^c} K_{beta,Omega}:
    cross_factor = ivq(4) if mutant == "p1-1-action-cross-factor" else ivq(2)
    phi_single_direct = cross_factor * pref * rect_int(ivq(0), delta, delta, beta)

    # Single-interval closed form: (2 U / Omega) sinh(Omega Delta / 2) sinh(Omega (beta - Delta) / 2) / sinh(beta Omega / 2)
    phi_single = (2 * u_val / omega) * iv_sinh(omega * delta / 2) * iv_sinh(
        omega * (beta - delta) / 2) / iv_sinh(beta * omega / 2)

    # Self-interaction integral 2 \iint_{[0,Delta]^2} K(s - u) ds du:
    self_int = 2 * pref * sq_int(delta)
    decomp_rhs = (u_val * delta - self_int) if mutant != "retardation-kernel-conservation" else (2 * u_val * delta)
    require(enclosed_zero(phi_single_direct - phi_single), "p1-1-action-cross-factor")
    require(enclosed_zero(phi_single - decomp_rhs), "retardation-kernel-conservation")
    require(bool(phi_single.a >= 0 and phi_single.b <= (u_val * delta).a), "retardation-kernel-conservation")


    # Finding R4: Verify the exact second derivative Phi_beta''(Delta) = -U * Omega * cosh(Omega*(beta/2 - Delta)) / sinh(beta*Omega/2)
    # (no spurious 1/2 factor) and the exact symmetric second-difference identity:
    def phi_of_d(d_arg):
        return (u_val / (omega * iv_sinh(beta * omega / 2))) * (
            iv_cosh(beta * omega / 2) - iv_cosh(omega * (beta / 2 - d_arg))
        )

    require(enclosed_zero(phi_of_d(delta) - phi_single), "r4-action-second-derivative")
    d2_factor = ivq(Q(1, 2)) if mutant == "r4-action-second-derivative" else ivq(1)
    phi_second_deriv = -d2_factor * u_val * omega * iv_cosh(omega * (beta / 2 - delta)) / iv_sinh(beta * omega / 2)
    h_step = ivq(Q(1, 1000))
    second_diff = (phi_of_d(delta + h_step) - 2 * phi_of_d(delta) + phi_of_d(delta - h_step)) / (h_step ** 2)
    sinc_sq = (iv_sinh(omega * h_step / 2) / (omega * h_step / 2)) ** 2
    require(
        enclosed_zero(second_diff - sinc_sq * phi_second_deriv)
        and bool(
            ivq(Q(-13116657148, 10 ** 10)).b < phi_second_deriv.a
            and phi_second_deriv.b < ivq(Q(-13116657147, 10 ** 10)).a
        ),
        "r4-action-second-derivative",
    )
    # Variance sandwich with corrected factor (beta * U / 2) * sum_x Var(n_x):
    var_sum = 2 * (delta / beta) * (1 - delta / beta)
    bh = beta * omega / 2
    var_pref = (beta * u_val) if mutant == "p1-1-action-cross-factor" else (beta * u_val / 2)
    var_low = (bh / iv_sinh(bh)) * var_pref * var_sum
    var_high = (bh * iv_cosh(bh) / iv_sinh(bh)) * var_pref * var_sum
    require(bool(var_low.a <= phi_single.b and phi_single.a <= var_high.b), "p1-1-action-cross-factor")

    # Two-interval two-state excursion A2 = [0, d1] U [s2, s2 + d2]:
    d1 = ivq(Q(2, 5))
    s2 = ivq(Q(2))
    d2 = ivq(Q(1, 2))
    meas_a2 = d1 + d2
    phi_two_direct = 2 * pref * (
        rect_int(ivq(0), d1, d1, s2)
        + rect_int(ivq(0), d1, s2 + d2, beta)
        + rect_int(d1, s2, s2, s2 + d2)
        + rect_int(s2, s2 + d2, s2 + d2, beta)
    )
    self_a2 = 2 * pref * (sq_int(d1) + sq_int(d2) + 2 * rect_int(ivq(0), d1, s2, s2 + d2))
    phi_two_decomp = u_val * meas_a2 - self_a2
    require(enclosed_zero(phi_two_direct - phi_two_decomp), "retardation-kernel-conservation")
    require(bool(phi_two_direct.a >= 0 and phi_two_direct.b <= (u_val * meas_a2).a),
            "retardation-kernel-conservation")

    # Verify the 3-state 4-hop closed walk DE -> OO -> ED -> OO -> DE counterexample (beta=4, U=1, |A|=1/5,
    # DE measure 19/10, ED measure 19/10): sum_x Var(n_x) = 19/10, so slow-limit Phi = (beta U / 2) * (19/10) = 19/5 > 2U|A| = 2/5 > U|A| = 1/5
    meas_a_3state = Q(1, 5)
    var_sum_3state = Q(19, 10)
    phi_slow_3state = Q(4 * 1, 2) * var_sum_3state
    require(phi_slow_3state == Q(19, 5) and phi_slow_3state > 2 * meas_a_3state > meas_a_3state,
            "retardation-kernel-conservation")

    # 2. Finding P1-2: Verify the reviewer's exact rational counterexample at U=1, beta=100, Omega=10^(-6), t=10^(-6),
    # centered mu=0 (z=1), proving that full physical w({e}) < 54e-12 while the erroneous doubled (4.16) gives > 65e-12,
    # whereas the corrected 4-family 2-hop lower bound (coefficient 1) gives ~33.22e-12 < w({e})!
    e10_low = sum(Q(10 ** k, factorial(k)) for k in range(40))
    require(e10_low > Q(20000), "p1-2-two-hop-double-count")
    x_bh = Q(1, 20000)
    require(1 - x_bh * x_bh > Q(99, 100), "p1-2-two-hop-double-count")
    t_ce = Q(1, 10 ** 6)
    beta_ce = Q(100)
    w2_coeff_upper = beta_ce * (Q(1000, 1891) + Q(1, 20000) + 100 * Q(1, 1 << 50))
    norm_x = 2 * beta_ce * t_ce
    tail_upper = (norm_x ** 4) / (24 * (1 - norm_x))
    full_w_upper = (t_ce ** 2) * w2_coeff_upper + tail_upper

    family_coeff = Q(2) if mutant == "p1-2-two-hop-double-count" else Q(1)
    # Using exp(-50) < 1/1000, 3*beta*U - 1 + exp(-300) > 299, and 2*beta*U - 1 + exp(-200) > 199:
    two_hop_lower_3u_ce = family_coeff * (t_ce ** 2) * Q(299) / (9 * ((1 + Q(1, 1000)) ** 2))
    two_hop_lower_2u_ce = family_coeff * (t_ce ** 2) * Q(199) / (4 * ((1 + Q(1, 1000)) ** 2))
    wrong_doubled_lower_ce = Q(2) * (t_ce ** 2) * Q(299) / (9 * ((1 + Q(1, 1000)) ** 2))
    require(full_w_upper < Q(54, 10 ** 12) < Q(65, 10 ** 12) < wrong_doubled_lower_ce,
            "p1-2-two-hop-double-count")
    require(two_hop_lower_3u_ce < two_hop_lower_2u_ce < full_w_upper, "p1-2-two-hop-double-count")

    # 3. Verify All-Order Single-Bond Lower Bound and Raw Mayer Expansion Obstruction Witness:
    # At U = 1, t = 1/4, beta = 640, centered mu = 0 (z = 1):
    # Corrected 4-family 2-hop bounds (kappa = 3 and kappa = 2) and all-order H_{2,kappa} bounds:
    u_nogo = ivq(Q(1))
    t_nogo = ivq(Q(1, 4))
    beta_nogo = ivq(Q(640))
    w_m1_simple = beta_nogo * (t_nogo ** 2) / (36 * u_nogo)
    a_nogo = mp.iv.exp(-beta_nogo * u_nogo / 2)
    denom = (1 + a_nogo) ** 2
    w_m1_exact = (t_nogo ** 2) * (3 * beta_nogo * u_nogo - 1 + mp.iv.exp(-3 * beta_nogo * u_nogo)) / (
        9 * (u_nogo ** 2) * denom
    )
    w_m1_2u_exact = (t_nogo ** 2) * (2 * beta_nogo * u_nogo - 1 + mp.iv.exp(-2 * beta_nogo * u_nogo)) / (
        4 * (u_nogo ** 2) * denom
    )
    delta_e3 = mp.iv.sqrt(9 * (u_nogo ** 2) + 4 * (t_nogo ** 2))
    lam_minus3 = (3 * u_nogo - delta_e3) / 2
    lam_plus3 = (3 * u_nogo + delta_e3) / 2
    c_minus3 = (3 * u_nogo + delta_e3) / (2 * delta_e3)
    c_plus3 = (delta_e3 - 3 * u_nogo) / (2 * delta_e3)
    h2_3u_00 = c_minus3 * mp.iv.exp(-beta_nogo * lam_minus3) + c_plus3 * mp.iv.exp(-beta_nogo * lam_plus3)
    w_allorder_low = (h2_3u_00 - 1) / denom

    delta_e2 = mp.iv.sqrt(4 * (u_nogo ** 2) + 4 * (t_nogo ** 2))
    lam_minus2 = (2 * u_nogo - delta_e2) / 2
    lam_plus2 = (2 * u_nogo + delta_e2) / 2
    c_minus2 = (2 * u_nogo + delta_e2) / (2 * delta_e2)
    c_plus2 = (delta_e2 - 2 * u_nogo) / (2 * delta_e2)
    h2_2u_00 = c_minus2 * mp.iv.exp(-beta_nogo * lam_minus2) + c_plus2 * mp.iv.exp(-beta_nogo * lam_plus2)
    w_allorder_2u_low = (h2_2u_00 - 1) / denom
    if mutant == "all-order-nogo-lower-bound":
        w_allorder_low = ivq(Q(1, 2))
    require(
        bool(
            w_m1_simple.a > 1
            and w_m1_2u_exact.a > w_m1_exact.b > 2 * w_m1_simple.b
            and w_allorder_2u_low.a > w_allorder_low.b > 500000
        ),
        "all-order-nogo-lower-bound",
    )

    # 4. Finding R6: Verify the exact 2x2 growth rate sqrt(U^2 + t^2) - U (and prefactor (1 + U/sqrt(U^2+t^2))/2)
    # for H_{2,2} at t = U = 1, where sqrt(2) - 1 < 1/2 = t^2 / (2U):
    u_r6, t_r6, beta_r6 = ivq(1), ivq(1), ivq(100)
    exact_rate_r6 = (
        (t_r6 ** 2) / (2 * u_r6)
        if mutant == "r6-exact-two-state-rate"
        else (mp.iv.sqrt(u_r6 ** 2 + t_r6 ** 2) - u_r6)
    )
    c_minus_r6 = (ivq(1) + u_r6 / mp.iv.sqrt(u_r6 ** 2 + t_r6 ** 2)) / 2
    c_plus_r6 = (ivq(1) - u_r6 / mp.iv.sqrt(u_r6 ** 2 + t_r6 ** 2)) / 2
    h22_100 = c_minus_r6 * mp.iv.exp(beta_r6 * (mp.iv.sqrt(ivq(2)) - ivq(1))) + c_plus_r6 * mp.iv.exp(
        -beta_r6 * (mp.iv.sqrt(ivq(2)) + ivq(1))
    )
    require(
        bool(
            ivq(Q(4142135, 10 ** 7)).b < exact_rate_r6.a
            and exact_rate_r6.b < ivq(Q(4142136, 10 ** 7)).a
            and (h22_100 / mp.iv.exp(ivq(50))).b < ivq(Q(1, 1000)).a
        ),
        "r6-exact-two-state-rate",
    )

    return phi_single, w_m1_exact, w_allorder_low


def certificate(mutant=None, verbose=True):
    require(mutant is None or mutant in MUTANTS, "known-mutant")
    n_probes = check_unconditional_sos(mutant)
    c7, q7, a7, q6 = check_fugacity_and_constants(mutant)
    theta7, theta4 = check_regime_map(mutant)
    phi_single, w_m1, w_all = check_retardation_and_nogo(mutant)

    if verbose:
        print(f"SOS_PROBES: {n_probes} exact 16-state rational probes across 5 (g,alpha) pairs up to g/alpha=1000")
        print("FUGACITY_CAP: 5*z_a - 4*z_b == (1-z)^2 + 2z(2w-1)^2 + 2zw^2 > 0; centered mu=0 (z=1) half-filling verified")
        print(f"UNCONDITIONAL_CONSTANTS: L=1024, C7={c7}=2^27, q7={q7}, A7={a7}, q_src={q6} (scalar source bound)")
        print(f"REGIME_CROSSOVERS: numerators 2^(89/2) and 2^28, asymptotic (2/3)^(1/3); at x=2^-19, theta_7={theta7} < theta_4={theta4}")
        print(f"RETARDATION_BARRIER: two-state Phi_beta in [{phi_single.a}, {phi_single.b}] <= U*Delta; "
              f"single-bond w({{e}}) >= {w_all.a} > {w_m1.a} > 1 uniformly in Omega>0")
        print("SCOPE: finite exact/interval witnesses; analytical theorems in RIGOR.md require independent review")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutant", choices=MUTANTS)
    parser.add_argument("--all-mutants", action="store_true")
    args = parser.parse_args()

    if args.mutant:
        try:
            certificate(mutant=args.mutant, verbose=False)
        except Rejected as exc:
            print(f"REJECTED {args.mutant}: {exc}")
            return 1
        print(f"UNCAUGHT {args.mutant}")
        return 0

    try:
        certificate(verbose=True)
    except Rejected as exc:
        print(f"BASELINE FAILURE: {exc}")
        return 1

    if not args.all_mutants:
        return 0

    script = str(Path(__file__).resolve())
    for mutant in MUTANTS:
        child = subprocess.run(
            [sys.executable, "-B", script, "--mutant", mutant],
            capture_output=True,
            text=True,
            check=False,
        )
        print(child.stdout.strip())
        if child.stderr:
            print(child.stderr.strip())
        print(f"mutant {mutant}: exit {child.returncode}")
        if child.returncode != 1 or not child.stdout.startswith(f"REJECTED {mutant}:"):
            return 1
    print(f"CERTIFIED: all {len(MUTANTS)} mutants rejected with exit 1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
