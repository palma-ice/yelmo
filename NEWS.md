# Yelmo release/tag notes

## Unreleased

- New `libs/fast_hydrology/` module providing distributed-water-flux / effective-pressure hydrology with method dispatch in the style of `basal_hydro_simple`. Currently exposes the K24 model (`method = 2`); `method = 0` is no-op (state frozen) and `method = 1` is reserved.
- K24 physics (constants, potential filling, smoothed gradients, flux accumulation, cavity-height blend, effective pressure) now lives in `fast_hydrology_k24` (`libs/fast_hydrology/k24.f90`). Constants are loaded at runtime from the `&fast_hydrology_k24` namelist group; defaults reproduce the previous hardcoded values.
- Two interchangeable flux solvers selectable via `flux_solver`: `0` = recursive depth-first (default; preserves current numerics), `1` = topological-sort (iterative, stack-safe; matches CISM's approach). FastHydrology.jl observed small numerical differences between the two — recursive remains default until the toposort variant is validated against a real run.
- API contract: the persistent water-layer thickness (`hyd%now%H_w`) is externally mutable by the host model; K24 reads ice geometry, bed, melt, sliding speed and Glen rate factor and returns `q_x`, `q_y`, `N`, `p_w` as diagnostics without touching `H_w`.
- New build dependency: FFTW3 (`-lfftw3`). Available via `fesm-utils`; configs that build this module need the link flag added.

## v1.15 (2026-01-19)

- Use of Gaussian Quadrature module from fesm-utils for calculating the Jacobian of velocity (strain-rate tensor), vertical velocity, dynamic viscosity (DIVA, SSA), basal friction, and other quantities.
- Added `uz_lim` to vertical velocity (improves stability for edge cases).
- Added new parameter `ytopo.dHdt_dyn_lim` to be able to limit rate of change of ice thickness due to dynamics can be (c
an help with stability).
- Implementation of switches to test different staggering methods (simple staggering versus Gaussian Quadrature, etc.): `ytherm.qb_method` and `ydyn.uz_method`.
- Implementation of LSF for calving, including CalvMIP test cases.
- Separation of calving parameters from topo parameters in namelist groups.
- Converted all further instances of get_neighbor_indices to get_neighbor_indices_bc_codes. Overall change led to spee
dup of 10% on a 16km Greeland run.
- OpenMP improvements means significant speedups are now possible for high-resolution runs.
- calc_bmb_total: bug fix; removed all traces of grounded_melt parameter, which was no longer used, and also removed optional argument mask_pd.
- Introduced `ydyn.scale_T` and `ydyn.T_frz` to control a linear reduction in friction until cf_ref in the case that ice is frozen at the base. This should make basal velocities more consistent with expectations, even when background friction is artificially low.
