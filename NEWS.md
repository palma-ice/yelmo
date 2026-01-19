# Yelmo release/tag notes

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
