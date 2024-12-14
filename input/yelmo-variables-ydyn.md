# ydyn

| id | variable          | dimensions       | units       | long_name                                     |
|----|-------------------|------------------|-------------|-----------------------------------------------|
|  1 | ux                | xc, yc, zeta     | m/yr        | x-velocity                                    |
|  2 | uy                | xc, yc, zeta     | m/yr        | y-velocity                                    |
|  3 | uxy               | xc, yc, zeta     | m/yr        | Horizonal velocity magnitude                  |
|  4 | uz                | xc, yc, zeta_ac  | m/yr        | z-component velocity                          |
|  5 | uz_star           | xc, yc, zeta_ac  | m/yr        | z-velocity with corr. for thermal advection   |
|  6 | ux_bar            | xc, yc           | m/yr        | Depth-averaged x-velocity                     |
|  7 | uy_bar            | xc, yc           | m/yr        | Depth-averaged y-velocity                     |
|  8 | uxy_bar           | xc, yc           | m/yr        | Depth-averaged horizontal velocity magnitude  |
|  9 | ux_bar_prev       | xc, yc           | m/yr        | Previous depth-averaged x-velocity            |
| 10 | uy_bar_prev       | xc, yc           | m/yr        | Previous depth-averaged y-velocity            |
| 11 | ux_b              | xc, yc           | m/yr        | Basal x-velocity                              |
| 12 | uy_b              | xc, yc           | m/yr        | Basal y-velocity                              |
| 13 | uz_b              | xc, yc           | m/yr        | Basal z-velocity                              |
| 14 | uxy_b             | xc, yc           | m/yr        | Basal horizontal velocity magnitude           |
| 15 | ux_s              | xc, yc           | m/yr        | Surface x-velocity                            |
| 16 | uy_s              | xc, yc           | m/yr        | Surface y-velocity                            |
| 17 | uz_s              | xc, yc           | m/yr        | Surface z-velocity                            |
| 18 | uxy_s             | xc, yc           | m/yr        | Surface horizontal velocity magnitude         |
| 19 | ux_i              | xc, yc, zeta     | m/yr        | Shearing x-velocity                           |
| 20 | uy_i              | xc, yc, zeta     | m/yr        | Shearing y-velocity                           |
| 21 | ux_i_bar          | xc, yc           | m/yr        | Depth-averaged shearing x-velocity            |
| 22 | uy_i_bar          | xc, yc           | m/yr        | Depth-averaged shearing y-velocity            |
| 23 | uxy_i_bar         | xc, yc           | m/yr        | Depth-averaged horizontal velocity magnitude  |
| 24 | duxydt            | xc, yc           | m/yr^2      | Time derivative of uxy                        |
| 25 | duxdz             | xc, yc, zeta     | 1/yr        | x-velocity vertical gradient                  |
| 26 | duydz             | xc, yc, zeta     | 1/yr        | y-velocity vertical gradient                  |
| 27 | duxdz_bar         | xc, yc           | 1/yr        | Depth-averaged x-velocity vertical gradient   |
| 28 | duydz_bar         | xc, yc           | 1/yr        | Depth-averaged y-velocity vertical gradient   |
| 29 | taud_acx          | xc, yc           | Pa          | Driving stress (x-dir)                        |
| 30 | taud_acy          | xc, yc           | Pa          | Driving stress (y-dir)                        |
| 31 | taud              | xc, yc           | Pa          | Driving stress magnitude                      |
| 32 | taub_acx          | xc, yc           | Pa          | Basal stress (x-dir)                          |
| 33 | taub_acy          | xc, yc           | Pa          | Basal stress (y-dir)                          |
| 34 | taub              | xc, yc           | Pa          | Basal stress magnitude                        |
| 35 | taul_int_acx      | xc, yc           | Pa          | Depth-integrated lateral stress (x-dir)                       |
| 36 | taul_int_acy      | xc, yc           | Pa          | Depth-integrated lateral stress (y-dir)                       |
| 37 | qq_gl_acx         | xc, yc           | m^3/yr      | Flux across grounding line                    |
| 38 | qq_gl_acy         | xc, yc           | m^3/yr      | Flux across grounding line                    |
| 39 | qq_acx            | xc, yc           | m^3/yr      | Flux (x-dir)                                  |
| 40 | qq_acy            | xc, yc           | m^3/yr      | Flux (y-dir)                                  |
| 41 | qq                | xc, yc           | m^3/yr      | Flux magnitude                                |
| 42 | de_eff            | xc, yc, zeta     | 1/yr        | Effective strain rate                         |
| 43 | visc_eff          | xc, yc, zeta     | Pa yr       | Effective viscosity                           |
| 44 | visc_eff_int      | xc, yc           | Pa yr m     | Depth-integrated viscosity                    |
| 45 | N_eff             | xc, yc           | Pa          | Effective pressure                            |
| 46 | cb_tgt            | xc, yc           | Pa          | Target basal parameter                        |
| 47 | cb_ref            | xc, yc           | --          | Reference basal parameter                     |
| 48 | c_bed             | xc, yc           | Pa          | Basal drag coefficient                        |
| 49 | beta_acx          | xc, yc           | Pa yr m^-1  | Basal stress factor (x)                       |
| 50 | beta_acy          | xc, yc           | Pa yr m^-1  | Basal stress factor (y)                       |
| 51 | beta              | xc, yc           | Pa yr m^-1  | Basal stress factor mag.                      |
| 52 | beta_eff          | xc, yc           | Pa yr m^-1  | Effective basal factor                        |
| 53 | f_vbvs            | xc, yc           | -           | Vertical basal stress                         |
| 54 | ssa_mask_acx      | xc, yc           | -           | SSA mask (x-dir)                              |
| 55 | ssa_mask_acy      | xc, yc           | -           | SSA mask (y-dir)                              |
| 56 | ssa_err_acx       | xc, yc           | m/yr        | SSA error (x-dir)                             |
| 57 | ssa_err_acy       | xc, yc           | m/yr        | SSA error (y-dir)                             |
| 58 | jvel_dxx          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duxdx             |
| 59 | jvel_dxy          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duxdy             |
| 60 | jvel_dxz          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duxdz             |
| 61 | jvel_dyx          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duydx             |
| 62 | jvel_dyy          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duydy             |
| 63 | jvel_dyz          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duydz             |
| 64 | jvel_dzx          | xc, yc, zeta_ac  | 1/yr        | Velocity Jacobian component duzdx             |
| 65 | jvel_dzy          | xc, yc, zeta_ac  | 1/yr        | Velocity Jacobian component duzdy             |
| 66 | jvel_dzz          | xc, yc, zeta_ac  | 1/yr        | Velocity Jacobian component duzdz             |
