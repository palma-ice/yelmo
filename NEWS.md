# Yelmo release/tag notes

## Unreleased

- New `libs/fast_hydrology/` module providing distributed-water-flux / effective-pressure hydrology with method dispatch in the style of `basal_hydro_simple`. Currently exposes the K24 model (`method = 2`); `method = 0` is no-op (state frozen) and `method = 1` is reserved.
- K24 physics (constants, potential filling, smoothed gradients, flux accumulation, cavity-height blend, effective pressure) now lives in `fast_hydrology_k24` (`libs/fast_hydrology/k24.f90`). Constants are loaded at runtime from the `&fast_hydrology_k24` namelist group; defaults reproduce the previous hardcoded values.
- Two interchangeable flux solvers selectable via `flux_solver`: `0` = recursive depth-first (default; preserves current numerics), `1` = topological-sort (iterative, stack-safe; matches CISM's approach). FastHydrology.jl observed small numerical differences between the two — recursive remains default until the toposort variant is validated against a real run.
- API contract: the persistent water-layer thickness (`hyd%now%H_w`) is externally mutable by the host model; K24 reads ice geometry, bed, melt, sliding speed and Glen rate factor and returns `q_x`, `q_y`, `N`, `p_w` as diagnostics without touching `H_w`.
- New build dependency: FFTW3 (`-lfftw3`). Available via `fesm-utils`; configs that build this module need the link flag added.
- New SSA momentum-balance formulation: **energy-based assembler** (`src/physics/solver_ssa_ac_energy.f90`), selectable via the new `ydyn.ssa_solver` parameter (`"residual"` | `"energy"`). The energy form assembles the symmetric positive-definite Hessian of the discrete viscous-energy density (membrane + shear + basal drag + driving), with η, β, H frozen per Picard iteration. Algebraic identity `K_inner = -A_residual_inner * dx*dy` for interior cells; calving-front BC encoded as a linear boundary-work term in the RHS rather than a one-sided FD stencil. SSA and DIVA both dispatch to the new assembler when `ssa_solver = "energy"`. The `ssa_lis_opt` namelist parameter has been split into `ssa_lis_opt_residual` and `ssa_lis_opt_energy` so each formulation carries its own LIS solver settings (defaults: `bicgsafe + jacobi` for residual, `cg + jacobi` for energy — CG is safe because the energy Hessian is SPD by construction).
- **initmip-grl finding (GRL-16KM, 100 yr, AB-SAM, dt_method=0)**: with `ssa_lat_bc="all"` (lateral BC applied to all ice fronts including grounded marine), the residual solver takes 48.9 s, the energy solver 5.8 s — **8.4× speedup**. Picard iterations: 2030 (residual) vs 135 (energy). The residual solver is BC-sensitive — switching from `"floating"` to `"all"` raises `uxy_s_max` from 4948 → 5424 m/yr and triples the Picard work. The energy solver is essentially insensitive to this switch (`uxy_s_max` literally identical, 3567 m/yr in both). This supports the prior intuition that the residual one-sided FD discretization of the calving-front BC is deficient at grounded marine fronts; the energy formulation's variational BC handling avoids the instability.
- New defaults: `ssa_solver = "energy"` in all parameter files; `ssa_lat_bc = "all"` in production files (`yelmo_initmip.nml`, `yelmo_TROUGH-F17.nml`, `yelmo_MISMIP+.nml`, `yelmo_calvingmip.nml`). Idealised infinite-slab benchmarks (`yelmo_SLAB-S06.nml`, `yelmo_ISMIPHOM.nml`) keep `ssa_lat_bc = "slab"` since their physics requires the slab extension.
- Algebraic-identity validation test: `tests/test_ssa_energy.f90` (build target `make ssa_energy`). Confirms the energy assembler matches `-A_residual * dx*dy` to single-precision roundoff on a periodic synthetic grid.

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
