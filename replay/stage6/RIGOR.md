# STAGE6: a uniform-frequency strong-coupling Holstein clustering theorem

## 1. Verdict, scope and explicit constants

**PROVED (analytical; independent review pending):** for the spinful square-lattice Holstein model, put `U=F²/K>0`, `Ω>0`, `τ=k_B T>0`, `β=1/τ`, and `t≥0`. There is no direct electron interaction. Define

```
L = 2^15 = 32768,
eps0 = 1/(64L) = 2^(-21),
C = 256L² = 2^38 = 274877906944,
q4 = 27/2048,
A4 = 9/(1-q4) = 18432/2021.
```

On square tori of side at least four, at every real chemical potential, if

```
t/U ≤ eps0,                  βt²/U ≤ 1/C,                 (1.1)
```

then, for distinct sites and graph distance `d=d(x,y)`, the original, equal-time electronic pair field `Δ_x=c_(x,down)c_(x,up)` satisfies

```
|<Δ_x† Δ_y>_GC| ≤ A4 q4^d.                                (1.2)
```

This is uniform in `Ω>0`, volume and chemical potential. Its constants are deliberately conservative, not optimized. In particular the on-site electronic pair-correlation matrix has uniformly summable rows; no nonexponential on-site pair order is possible. The same exponential bound holds in translation-invariant thermodynamic limits of total-number canonical states at every fixed interior density, and in balanced separately fixed spin-number sectors, by the source argument in §7. A uniform finite-volume canonical estimate and imbalanced spin-sector estimates are **not** claimed.

Consequently, with pair order defined by these on-site electronic correlations,

```
τ_pair ≤ C t²/U,                 t/U ≤ eps0.              (1.3)
```

No analogous bond-SSH theorem is established: §9 explains the missing structural prerequisite. No claim is made about the full unrestricted Yang matrix, experimental transition temperatures, or a linear-in-Ω ceiling. Combining (1.3) with the released LINEAR4 Holstein bound gives the minimum of their two thresholds. This can strictly improve the previous upper bound at finite parameters (§8), but does not improve LINEAR4's fixed-(t,U) adiabatic exponent.

**CERTIFIED** in the replay means finite exact rational or outward interval arithmetic only; it is not machine verification of this analytical many-body theorem. No external published theorem is used as authority.

## 2. Exact chronological contours and spatial support factorization

Use the centered Holstein coupling and its on-site Lang-Firsov unitary. Subtract the oscillator zero-point constant. The transformed atomic Hamiltonian is

```
H0 = Σ_x [Ω a_x†a_x -(U/2)(n_x-1)² - μ n_x].
```

An additive `-U/2` per site cancels from every normalized trace. Thus its electronic part can instead be written as

```
D0 = (U/2) Nodd - μ N,
Nodd = Σ_x 1_(n_x=1),           N = Σ_x n_x.               (2.1)
```

The displacement per electron hop has `α²=U/(2Ω)`. Each directed spin hop is `t` times a fermionic partial isometry and endpoint displacements. The full transformed bond hopping has norm at most `2t` and is fermion-even.

For a closed Fock trajectory, let `n_x(u)` be its periodic, piecewise constant local occupation, including jumps caused by any inserted physical pair source. Its exact normalized oscillator trace is

```
Tr_ph,0[ T prod_j D_x(c_j α,u_j) ] = exp(-Φβ[n_x]),
T prod_j D_x(c_j α,u_j) = D_x(c_m α,u_m)...D_x(c_1 α,u_1),
0≤u_1<...<u_m<β,                Σ_j c_j=0,
Φβ[n] = (1/4) ∫_0^β du ∫_0^β dv (n(u)-n(v))² Kβ(u-v),
Kβ(w) = (UΩ/2) [e^(-Ωw)+e^(-Ω(β-w))]/[1-e^(-βΩ)],
                                      0≤w≤β.             (2.2)
```

The kernel is periodically extended. Therefore `Kβ>0`, `∫_0^β Kβ=U`, and `Φβ≥0`. In particular the complete phonon factor, with the physical pair insertions included, belongs to `(0,1]`.

For clarity, the oscillator calculation is not an increasing-time product. If `D(cα,u)=exp(cα e^(Ωu)a†-cα e^(-Ωu)a)`, scalar BCH and thermal normal ordering give

```
log Tr_ph,0(T prod D)
 = -(α²/2)coth(βΩ/2)(Σc_j)²
   - α² Σ_(i<j) c_i c_j
       [cosh(βΩ/2-Ω(u_j-u_i))-cosh(βΩ/2)]/sinh(βΩ/2).
```

For neutral jumps this equals `-Φβ`: integrate the exponential kernel twice by parts across the jumps of the periodic occupation. Positivity also follows immediately from the square in (2.2). These are the independently corrected identities in `LINEAR5-REVIEW/REVIEW-STAGE5.md` §2. Nothing here uses the failed STAGE5 pair envelope, failed channel norms, or excited-spectrum SW interpretation.

Dyson terms are grouped by their distinct used bond labels `E`. Connected components with disjoint endpoint sets act on independent oscillators and independent atomic electronic spaces. All bond and pair-source operators are even, so they commute on disjoint sets. Summing chronological shuffles gives an exact hard-core spatial polymer gas, whose polymers are nonempty connected edge-label sets and whose compatibility is disjoint vertex support. This is the same representation as the reviewed STAGE4 §4.1.

**Important:** no factorization in time is asserted. Every oscillator cross-action and every shared-site occupation correlation within a polymer is retained before any bound. Thus the static-site-variance problem is absorbed into the full spatial support rather than being dropped or conditioned away.

## 3. A signless finite-dimensional trajectory majorant

Fix a finite physical edge set `E` with endpoint set `S`, `v=|S|`, and `k=|E|`. Use the occupation basis with four states per site: empty, up, down, double. Define `B_e` to have a matrix entry one for each legal directed spin hop on edge `e`, with the CAR sign removed. The destination is unique for a given directed spin hop and initial state; distinct such hops do not give the same destination. Thus `B=Σ_e B_e` is real symmetric and entrywise nonnegative. It is a bookkeeping matrix, not a claimed physical Hamiltonian.

Let `z_a=1+z²+2z e^(-βU/2)`, `z=e^(βμ)>0`, be the shifted atomic electronic partition per site. For a scale `l≥1`, define

```
P_E(l) = Tr_el exp[-β(D0-l t B)] / z_a^v.                 (3.1)
```

Its Dyson coefficients are nonnegative: diagonal atomic heat factors and nonnegative hopping matrices suffice, with no fermionic signs. For each specified closed trajectory the absolute physical CAR factor is exactly the product of its hopping amplitudes; its remaining oscillator factor is at most one by (2.2). Consequently the sum `Wabs(E)` of absolute physical trajectory weights using every edge in `E` obeys

```
|w(E)| ≤ Wabs(E) ≤ l^(-k) P_E(l).                         (3.2)
```

Indeed a term covering `k` distinct labels has at least `k` hopping insertions; scaling every hop by `l` multiplies it by `l^m≥l^k`, and dropping the covering restriction only enlarges the positive sum. This is the missing alternative to trying to assign an independent `t/U` factor to each chronological gap.

The same bound holds for the marked pair-polymer activity. Strip CAR signs from the electronic matrix of `Δ_x†Δ_y` to obtain a nonnegative partial isometry `O_xy` with norm one. It is included at time zero in every trajectory, so jumps of magnitude two and their oscillator displacements remain inside (2.2). Nonzero terms still have neutral jumps at each site. After dropping their phonon factors, their positive sum is `Tr[O_xy exp(-β(D0-ltB))]/z_a^v`. For any positive matrix `Q` and contraction `O`,

```
|Tr(OQ)| ≤ Tr Q,
```

by diagonalizing `Q` and bounding each diagonal matrix element of `O` by one. Thus (3.2) also majorizes this marked sum. This does **not** upper-bound a correlation by dividing two unrelated positive partition functions: the physical denominator will be handled by the exact polymer deletion ratio in §6.

## 4. Form estimate: absorbing odd-sector hopping into half the atomic gap

**PROVED lemma.** On any finite physical subgraph of maximum degree `D`, let `g≥0` and assume `g/U≤1/(16D)`. Then, as quadratic forms on the finite electronic space,

```
g B ≤ (U/4) Nodd + (16D k g²/U) I.                       (4.1)
```

Here `k` counts physical edges only. Chemical potential plays no role, since hopping conserves `N`.

### Proof with the endpoint and spin factors retained

Classify every legal hop by its change in `Nodd`: pair creation `+2`, pair annihilation `-2`, or propagation `0`. Pair creation connects an even-even edge state DE or ED to an odd-odd state. There are at most two outgoing creation hops per edge from any configuration, hence at most `2k` globally. Pair annihilation requires two odd endpoints; there are at most two such hops per odd-odd edge. Since each odd vertex has at most `D` neighbors, its total outgoing annihilation degree is at most `D Nodd`.

Orient each creation/annihilation matrix edge from the lower-odd configuration `a` to the higher-odd configuration `b`. For any real or complex vector `ψ`, and any `h>0`,

```
2g |ψ_a ψ_b| ≤ g h |ψ_a|² + (g/h)|ψ_b|².
```

Sum over these oriented edges. The two degree counts give

```
g(B+ + B-) ≤ 2k g h I + (Dg/h) Nodd.                     (4.2)
```

A zero-channel hop touches at least one odd endpoint and preserves `Nodd`. The outgoing degree is bounded by `2D Nodd` (a deliberately loose bound). The elementary inequality `2|ψ_a ψ_b|≤|ψ_a|²+|ψ_b|²` and the symmetric adjacency row sums then give

```
g B0 ≤ 2D g Nodd.                                       (4.3)
```

For `g>0`, choose `h=8Dg/U`. Then `Dg/h=U/8`, `2Dg≤U/8`, and `2kgh=16Dkg²/U`, proving (4.1). At `g=0` it is immediate by continuity. In particular we never assign the incorrect norm `sqrt(2)t` to an entire pair-breaking channel.

## 5. The frequency-independent support bound

By (4.1),

```
D0-gB ≥ (U/4) Nodd - μN -(16Dkg²/U) I.
```

Eigenvalue monotonicity of Hermitian matrices, followed by the decreasing scalar function `exp(-βx)` and summation, implies trace monotonicity. Operator monotonicity of the matrix exponential is **not** being assumed. Hence

```
P_E(l) ≤ exp(16D k βl²t²/U) (z_b/z_a)^v,
z_b = 1+z²+2z e^(-βU/4).                                (5.1)
```

Writing `b=βU`, the inequality `2z≤1+z²` gives, uniformly for all real chemical potentials,

```
z_b/z_a ≤ 1+e^(-b/4) ≤ exp(e^(-b/4)) ≤ e.                (5.2)
```

Combine (3.2)-(5.2), choose `D=4`, `l=L`, and impose (1.1). Then `Lt/U≤1/64=1/(16D)` and `16D βL²t²/U≤1/4`. For every connected physical `E`, `v≤k+1`, and for both unmarked and marked activities,

```
|w(E)|, |w_xy(E)| ≤ L^(-k) exp(k/4+v).                   (5.3)
```

This is a bound on full, overlapping, all-order spatial supports. There are no initial-even-state restrictions, discarded initial spinons, short-time restrictions, resonance denominators, phonon cutoffs, or hidden `βt` assumptions. Thermally winding odd-electron trajectories are included in the finite matrix trace (3.1).

## 6. Cluster convergence and original pair correlations

The rooted number of connected distinct-edge sets of size `k` on a graph of degree `Dg` is at most `Dg^(2k)` by the doubled-edge Euler-tour injection (STAGE4 §4.2). For physical degree four, (5.3) implies

```
Σ_(E:x∈S_E) |w(E)| e^|S_E|
 ≤ e² Σ_(k≥1) [16 e^(9/4)/L]^k
 < 9 Σ_(k≥1) q4^k
 = 243/2021 < 1,                                         (6.1)
```

because `e<3`, `e^(9/4)<e³<27`, and `q4=16·27/L=27/2048`. The same estimate holds for marked supports. The value `e<3` follows from its positive Taylor series and the geometric tail estimate, with a rational enclosure in the replay.

For completeness the gas criterion used is: if `Σ_(γ:x∈Sγ)|wγ|e^|Sγ|≤1` at every vertex, its log cluster series converges absolutely and its support-deletion ratio satisfies

```
|Z_(Λ\S)/Z_Λ| ≤ e^|S|.                                  (6.2)
```

Its proof is the hard-core alternating-graph/tree bound in reviewed STAGE4 §4.3: partition connected incompatibility graphs by their minimum spanning tree; the alternating sum over the optional larger chords has magnitude zero or one. Rooting the remaining trees, their child sum is bounded by the recurrence `F_next(γ)=|wγ| exp(Σ_(γ' incompatible γ)F(γ'))`; the site criterion inductively gives `F(γ)≤|wγ|e^|Sγ|`. Clusters touching `S` have total log magnitude at most `|S|`, and exponentiating their deletion proves (6.2). This argument allows signed physical activities.

In the physical marked trace, a connected component containing only one pair insertion has zero trace by charge conservation on its support. Hence one marked polymer connects `x` and `y`. It has at least `d(x,y)` physical edges. Unmarked components outside it supply precisely the physical deletion ratio, not a signless denominator. Using (5.3), (6.2) and the same rooted counting gives

```
|<Δ_x†Δ_y>_GC|
 ≤ Σ_(E connected:x,y∈S_E) |w_xy(E)| e^|S_E|
 ≤ 9 q4^d/(1-q4) = A4 q4^d.
```

This proves (1.2). On a square torus a distance shell has at most `4d` sites. The on-site diagonal pair expectation is at most one, and `Σ_(d≥1) d q4^d=q4/(1-q4)²`; row summability is uniform in volume.

## 7. Real-source transfer to fixed-density canonical limits

This section transfers the proved bound, not an unproved finite-volume equivalence of charged correlators. Fix a nonzero displacement `r` and phase `θ`, and add

```
H(ξ)=H-ξ X,
X=Σ_x [e^(iθ) Δ_x†Δ_(x+r)+e^(-iθ) Δ_(x+r)†Δ_x].
```

Each source label is supported on two sites and has norm at most two. Count source labels as a multigraph if two labels share the same endpoints; the graph degree is at most six. Physical hopping still has maximum degree four. Source pair jumps preserve local occupation parity. After stripping signs/phases, the source bookkeeping matrix on a label has norm at most two; consequently the form estimate (4.1) changes only by an additive `2L|ξ|` per source label. In the complex disk

```
|ξ|≤ξ0 = τ/(16L),                                       (7.1)
```

its scaled trace penalty is at most `2βL|ξ|=1/8` per source label. Every physical hopping label costs at most `1/4` in (5.1). Thus (5.3) remains valid with `k` the total number of physical/source labels and `v≤k+1`.

The degree-six site budget is bounded by

```
q6=36·27/L=243/8192,
9 q6/(1-q6)=2187/7949 < 1.                               (7.2)
```

It follows that the physical log series is analytic uniformly in volume and chemical potential in (7.1). Its per-site absolute cluster budget is at most one; the normalized log partition difference from zero source is at most two. Cauchy's bounds give a remainder at most `4ξ²/ξ0²` for real `|ξ|≤ξ0/2`. The zero-source derivative is `2β Re(e^(iθ)<Δ_0†Δ_r>_GC)`. Therefore

```
V^(-1)[log Z_GC(ξ;μ)-log Z_GC(0;μ)]
 ≤ 2β A4 q4^d |ξ|+4ξ²/ξ0².                              (7.3)
```

The number-conserving source permits exactly the thermodynamic pressure argument in reviewed STAGE4 §5.3. Its hypotheses still hold: all transformed bonds and sources are bounded; for each fixed `r` the boundary interaction count is subextensive; atomic sector traces are finite; and the electron occupations per site are finite. To spell out the ensemble step, divide a large torus into fixed open blocks, dropping boundary bonds at an `O(boundary)` log-trace cost. The disconnected canonical partition is a sum over block number allocations. Its log differs from the largest allocation log by at most the number of blocks times `log(2 block_volume+1)`. The limiting canonical entropy/free-energy density is the concave envelope of the finite block-sector values as block size tends to infinity. Boundary cost and allocation entropy per volume vanish. This yields the concave canonical density function `f(n)` and its grand-canonical dual `p(μ)=max_n[f(n)+βμn]`, along arbitrary `N/V→n` sequences. At every interior density choose a finite supporting chemical potential, so

```
δ_V = V^(-1)[log Z_GC(0;μ)-log Z_N(0)-βμN] ≥0,
δ_V →0.                                                 (7.4)
```

Since `e^(βμN)Z_N(ξ)≤Z_GC(ξ;μ)`, combine (7.3)-(7.4). Convexity of the finite canonical real-source log partition bounds its derivative between secants at `±ξ`, giving

```
|2β Re(e^(iθ)<Δ_0†Δ_r>_N)|
 ≤ 2β A4 q4^d +4ξ/ξ0²+δ_V/ξ,        0<ξ≤ξ0/2.
```

Take the thermodynamic limit first at fixed `ξ`, then `ξ→0`, and vary `θ`. This proves the same exponential bound for each fixed displacement in translation-invariant canonical limits. No uniform finite-volume canonical deletion ratio is inferred.

The balanced separately fixed spin sectors also follow, rather than being assumed equivalent to total-number sectors. Repeat the block allocation argument with the pair `(N_up,N_down)`; allocation entropy per block is at most `2 log(block_volume+1)`, still vanishing per site. It gives a finite concave density function `f(n_up,n_down)` on the density square and a two-variable pressure dual. At a balanced interior point `(a,a)`, choose a finite supporting vector. Spin-exchange symmetry makes its swapped vector supporting as well, so their average is supporting and has equal components. Thus a common, real chemical potential supports `(a,a)` in the unrestricted grand-canonical ensemble. Equation (7.4) holds with `Z_N` replaced by `Z_(N_up,N_down)` and `N=N_up+N_down`, along balanced limiting density sequences. The source preserves both spin numbers, so the same secant proof applies without modification. This does not provide supporting common chemical potentials at imbalanced densities.

## 8. What is genuinely improved

Let `τ4` denote the positive LINEAR4 threshold determined by

```
τ4³=2^25 UΩ² [1+log(1+E/τ4)],
```

where `E` is its released model-dependent kinetic/coupling scale. In the strong-coupling corner proved here,

```
τ_pair ≤ min(Ct²/U, τ4).                                 (8.1)
```

**CERTIFIED exact finite-parameter comparison:** choose `t/U=2^(-30)` and `Ω/U=1`. These satisfy `t/U≤eps0`. The new ceiling is `Ct²/U=2^(-22)U`. STAGE4's bare-hopping ceiling `1024t` is `2^(-20)U`, four times larger. LINEAR4 has `τ4³≥2^25 U³`, hence `τ4>U`, so the STAGE6 ceiling is strictly smaller without needing to estimate its logarithm. This is a strict original-model improvement in an extremely strong-coupling family; the constants make no practical sharpness claim.

At fixed positive `(t,U)` in this corner, `Ct²/U` remains positive while `τ4→0` as `Ω→0`. Thus this theorem neither improves the adiabatic power nor closes the original linear-in-phonon-energy target. It closes the uniform-frequency **strong-coupling** on-site Holstein question stated in (1.1), with the ensemble qualifications above.

## 9. SSH and remaining qualifications

The proof uses the exact local identity `D0=(U/2)Nodd-μN` and an on-site Lang-Firsov transform with independent site displacement factors. Bond SSH instead couples bond electronic operators to independent bond oscillators. Neighboring bond operators do not commute, so no on-site occupation-diagonal transform or the form estimate (4.1) has been established for its atomic reference. The Holstein signless majorant cannot be reused for SSH merely by relabeling sites as bonds. **OPEN:** a separate SSH theorem with comparable uniform-frequency constants.

**OPEN beyond this result:** linear/polylog-linear frequency ceiling, unrestricted pair-matrix order, finite-volume uniform canonical correlations, imbalanced spin-sector canonical transfer, optimized constants, and a theorem outside strong coupling. Independent adversarial review of STAGE6 is pending.

## 10. References to verify and replay scope

No external paper or published theorem is assumed. Internal mathematical inputs are:

1. `LINEAR5-REVIEW/REVIEW-STAGE5.md` §2: chronological oscillator trace, positive contour kernel and `∫K=U` (independently corrected review).
2. `ADIABATIC/LINEAR5/STAGE4/RIGOR.md` §§4.1-4.4 and 5.3: exact hard-core support grouping, connected-edge counting, tree/deletion criterion, and thermodynamic source transfer. Its independent review is `LINEAR5-REVIEW/REVIEW-STAGE4.md`.
3. Released `ADIABATIC/LINEAR4/RESULT.md`: threshold used only for comparison in §8.

The replay checks an independently assembled 16-state CAR bond, its signless channel degrees and form estimate by exact rational diagonal-dominance witnesses, the exact fugacity-ratio algebra, chronological versus reverse-order oscillator witnesses by outward interval arithmetic, all explicit activity/source budgets and the finite-parameter threshold comparison. Guarded mutants corrupt each logical ingredient and must exit one. Finite arithmetic does not certify the infinite oscillator trace identity, arbitrary-graph analytical form proof, cluster expansion, or canonical thermodynamic limits. No Lean file is produced because these analytical arguments have not been formalized; a token algebra file would not substantiate them.
