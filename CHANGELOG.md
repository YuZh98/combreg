# combreg 0.2.0

### Added
- `crr_control()` gains an adaptive latent-utility block size. `zeta_block` now accepts `"adaptive"` in addition to a fixed integer. 
- `crr()` accepts a one-sided formula and `data` in place of a covariate matrix. The matrix interface is unchanged.
- `crr_constraints()` gains a `dedup` argument (default `TRUE`) that removes redundant constraint rows. 

### Changed
- `accept_rate` now averages over the post-warmup sampling phase only, so it reflects the frozen proposal rather than warmup adaptation transients.


# combreg 0.1.0

Initial release.