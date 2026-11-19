# Replay supplement for Strong-coupling suppression of electronic pair order in the Holstein model, uniformly in phonon frequency and without a hopping cutoff

Davis Rens, SuperMatics Inc. ([ORCID](https://orcid.org/0000-0001-5839-5197)).

Paper DOI: [10.5281/zenodo.24380128](https://doi.org/10.5281/zenodo.24380128). Supplement DOI: [10.5281/zenodo.24380140](https://doi.org/10.5281/zenodo.24380140).

## Scope and evidence

For the original spinful square-lattice Holstein model with $U=F^2/K>0$ and $t\ge0$, the manuscript proves exponential decay of equal-time on-site electronic pair correlations whenever $\beta t^2/U\le2^{-27}$. No restriction on $t/U$ is imposed. For distinct sites,

$$
|\langle\Delta_x^\dagger\Delta_y\rangle|
\le\frac{240}{59}\left(\frac{5}{64}\right)^{d(x,y)}.
$$

The bound is uniform in positive phonon frequency, real common chemical potential, and periodic square tori of side at least four. Its pair-order consequence is $k_B T_{\rm pair}\le2^{27}t^2/U$. Canonical transfer concerns translation-invariant thermodynamic limits at interior total density and balanced interior separately fixed spin densities, not uniform finite-volume canonical estimates.

This is a strong-coupling pair-order bound, not a linear-in-phonon-frequency ceiling. No material prediction, synthesis, or measured superconducting transition is claimed.

| Result | Label | Evidence |
|---|---|---|
| Uniform-frequency pair-clustering theorem and qualified consequences | PROVED | Admitted manuscript argument and companion proof notes, independently reviewed by separate AI agents |
| Specified finite exact and interval witnesses | CERTIFIED | Frozen replay, default and all 27 guarded mutants checked in the release CI |
| Linear-in-phonon-frequency pairing ceiling | OPEN | Not proved by this paper |

The replay checks finite algebraic and interval witnesses. It does not formalize the oscillator trace identity, infinite-volume cluster expansion, or canonical pressure argument. No empirical data or simulation results are used.

## Review status

The proof notes were reviewed in `REVIEW-STAGE7R` and `REVIEW-STAGE7R-DELTA`; the folded manuscript received PASS in `REVIEW-STAGE6-PAPER-FOLD` with no required fixes. All these reviews were carried out by separate AI agents. Historical labels retained in the bundled notes predate the completed reviews. This preprint has not been peer reviewed.

## Contents and replay

The manuscript is in `paper/main.tex` and `paper/refs.bib`. The byte-exact finite replay and its original source pin are `replay/replay.py` and `replay/SOURCE-PINS.json`. The strengthened proof and its result note are `replay/stage7r/RIGOR.md` and `replay/stage7r/RESULT.md`. The predecessor `replay/stage6/RIGOR.md` is included only as a historical source cited by the strengthened note; its cutoff and constants are superseded, not the current theorem. `ADMITTED.sha256` pins every frozen input, including the pre-metadata manuscript and PDF hashes. `MANIFEST.sha256` covers the released tree.

```sh
python -m pip install mpmath==1.4.1
python -B replay/replay.py
python -B replay/replay.py --all-mutants
sha256sum -c MANIFEST.sha256
(cd paper && tectonic --keep-logs main.tex)
```

## Licenses and provenance

The replay code and documentation use Apache-2.0. The manuscript text in `paper/` is distributed under CC-BY-4.0 through its paper record. AI agents generated proofs, certificates, code, and manuscript drafts and edits under D. Rens' direction and review.
