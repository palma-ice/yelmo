# ydyn

| variable          | dimensions       | units       | long_name                                     |
|-------------------|------------------|-------------|-----------------------------------------------|
| ux                | xc, yc, zeta     | m/yr        | x-velocity                                    |
| uy                | xc, yc, zeta     | m/yr        | y-velocity                                    |
| uxy               | xc, yc, zeta     | m/yr        | Horizonal velocity magnitude                  |
| uz                | xc, yc, zeta_ac  | m/yr        | z-component velocity                          |
| uz_star           | xc, yc, zeta_ac  | m/yr        | z-velocity with corr. for thermal advection   |
| ux_bar            | xc, yc           | m/yr        | Depth-averaged x-velocity                     |
| uy_bar            | xc, yc           | m/yr        | Depth-averaged y-velocity                     |
| uxy_bar           | xc, yc           | m/yr        | Depth-averaged horizontal velocity magnitude  |
| ux_bar_prev       | xc, yc           | m/yr        | Previous depth-averaged x-velocity            |
| uy_bar_prev       | xc, yc           | m/yr        | Previous depth-averaged y-velocity            |
| ux_b              | xc, yc           | m/yr        | Basal x-velocity                              |
| uy_b              | xc, yc           | m/yr        | Basal y-velocity                              |
| uz_b              | xc, yc           | m/yr        | Basal z-velocity                              |
| uxy_b             | xc, yc           | m/yr        | Basal horizontal velocity magnitude           |
| ux_s              | xc, yc           | m/yr        | Surface x-velocity                            |
| uy_s              | xc, yc           | m/yr        | Surface y-velocity                            |
| uz_s              | xc, yc           | m/yr        | Surface z-velocity                            |
| uxy_s             | xc, yc           | m/yr        | Surface horizontal velocity magnitude         |
| ux_i              | xc, yc, zeta     | m/yr        | Shearing x-velocity                           |
| uy_i              | xc, yc, zeta     | m/yr        | Shearing y-velocity                           |
| ux_i_bar          | xc, yc           | m/yr        | Depth-averaged shearing x-velocity            |
| uy_i_bar          | xc, yc           | m/yr        | Depth-averaged shearing y-velocity            |
| uxy_i_bar         | xc, yc           | m/yr        | Depth-averaged horizontal velocity magnitude  |
| duxydt            | xc, yc           | m/yr^2      | Time derivative of uxy                        |
| duxdz             | xc, yc, zeta     | 1/yr        | x-velocity vertical gradient                  |
| duydz             | xc, yc, zeta     | 1/yr        | y-velocity vertical gradient                  |
| duxdz_bar         | xc, yc           | 1/yr        | Depth-averaged x-velocity vertical gradient   |
| duydz_bar         | xc, yc           | 1/yr        | Depth-averaged y-velocity vertical gradient   |
| taud_acx          | xc, yc           | Pa          | Driving stress (x-dir)                        |
| taud_acy          | xc, yc           | Pa          | Driving stress (y-dir)                        |
| taud              | xc, yc           | Pa          | Driving stress magnitude                      |
| taub_acx          | xc, yc           | Pa          | Basal stress (x-dir)                          |
| taub_acy          | xc, yc           | Pa          | Basal stress (y-dir)                          |
| taub              | xc, yc           | Pa          | Basal stress magnitude                        |
| taul_int_acx      | xc, yc           | Pa          | Depth-integrated lateral stress (x-dir)                       |
| taul_int_acy      | xc, yc           | Pa          | Depth-integrated lateral stress (y-dir)                       |
| qq_gl_acx         | xc, yc           | m^2/yr      | Flux across grounding line                    |
| qq_gl_acy         | xc, yc           | m^2/yr      | Flux across grounding line                    |
| qq_acx            | xc, yc           | m^2/yr      | Flux (x-dir)                                  |
| qq_acy            | xc, yc           | m^2/yr      | Flux (y-dir)                                  |
| qq                | xc, yc           | m^2/yr      | Flux magnitude                                |
| de_eff            | xc, yc, zeta     | 1/yr        | Effective strain rate                         |
| visc_eff          | xc, yc, zeta     | Pa yr       | Effective viscosity                           |
| visc_eff_int      | xc, yc           | Pa yr m     | Depth-integrated viscosity                    |
| N_eff             | xc, yc           | Pa          | Effective pressure                            |
| cb_tgt            | xc, yc           | Pa          | Target basal parameter                        |
| cb_ref            | xc, yc           | --          | Reference basal parameter                     |
| c_bed             | xc, yc           | Pa          | Basal drag coefficient                        |
| beta_acx          | xc, yc           | Pa a m^-1   | Basal stress factor (x)                       |
| beta_acy          | xc, yc           | Pa a m^-1   | Basal stress factor (y)                       |
| beta              | xc, yc           | Pa a m^-1   | Basal stress factor mag.                      |
| beta_eff          | xc, yc           | Pa a m^-1   | Effective basal factor                        |
| f_vbvs            | xc, yc           | -           | Vertical basal stress                         |
| ssa_mask_acx      | xc, yc           | -           | SSA mask (x-dir)                              |
| ssa_mask_acy      | xc, yc           | -           | SSA mask (y-dir)                              |
| ssa_err_acx       | xc, yc           | m/yr        | SSA error (x-dir)                             |
| ssa_err_acy       | xc, yc           | m/yr        | SSA error (y-dir)                             |
