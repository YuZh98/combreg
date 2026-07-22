# combreg (development version)

### Added
- `crr_control()` gains `n_iter_hit_and_run_mh`, setting the hit-and-run sweep count for the dual certificate drawn inside the Metropolis-Hastings update of the latent utilities, separately from the dual refresh. It defaults to `n_iter_hit_and_run`, so results are unchanged unless it is set.


# combreg 0.2.0

### Added
- `crr_control()` gains an adaptive latent-utility block size. `zeta_block` now accepts `"adaptive"` in addition to a fixed integer. 
- `crr()` accepts a one-sided formula and `data` in place of a covariate matrix. The matrix interface is unchanged.
- `crr_constraints()` gains a `dedup` argument (default `TRUE`) that removes redundant constraint rows. 

### Changed
- `accept_rate` now averages over the post-warmup sampling phase only, so it reflects the frozen proposal rather than warmup adaptation transients.


# combreg 0.1.0

Initial release.