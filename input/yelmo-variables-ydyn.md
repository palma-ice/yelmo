# ydyn

| variable          | dimensions       | units       | long_name                                     |
|-------------------|------------------|-------------|-----------------------------------------------|
| ux                | xc, yc, zeta     | m/s         | x-component velocity                          |
| uy                | xc, yc, zeta     | m/s         | y-component velocity                          |
| uxy               | xc, yc, zeta     | m/s         | xy-plane velocity                             |
| uz                | xc, yc, zeta     | m/s         | z-component velocity                          |
| uz_star           | xc, yc, zeta     | m/s         | Adjusted z-velocity                           |
| ux_bar            | xc, yc           | m/s         | Mean x-velocity                               |
| uy_bar            | xc, yc           | m/s         | Mean y-velocity                               |
| uxy_bar           | xc, yc           | m/s         | Mean xy-plane velocity                        |
| ux_bar_prev       | xc, yc           | m/s         | Previous mean x-velocity                      |
| uy_bar_prev       | xc, yc           | m/s         | Previous mean y-velocity                      |
| ux_b              | xc, yc           | m/s         | Basal x-velocity                              |
| uy_b              | xc, yc           | m/s         | Basal y-velocity                              |
| uz_b              | xc, yc           | m/s         | Basal z-velocity                              |
| uxy_b             | xc, yc           | m/s         | Basal xy-plane velocity                       |
| ux_s              | xc, yc           | m/s         | Surface x-velocity                            |
| uy_s              | xc, yc           | m/s         | Surface y-velocity                            |
| uz_s              | xc, yc           | m/s         | Surface z-velocity                            |
| uxy_s             | xc, yc           | m/s         | Surface xy-plane vel.                         |
| ux_i              | xc, yc, zeta     | m/s         | Interpolated x-velocity                       |
| uy_i              | xc, yc, zeta     | m/s         | Interpolated y-velocity                       |
| ux_i_bar          | xc, yc           | m/s         | Mean interpolated x-vel.                      |
| uy_i_bar          | xc, yc           | m/s         | Mean interpolated y-vel.                      |
| uxy_i_bar         | xc, yc           | m/s         | Mean interp. xy-velocity                      |
| duxydt            | xc, yc           | m/s^2       | Time derivative of uxy                        |
| duxdz             | xc, yc, zeta     | 1/s         | x-velocity vertical grad.                     |
| duydz             | xc, yc, zeta     | 1/s         | y-velocity vertical grad.                     |
| duxdz_bar         | xc, yc           | 1/s         | Mean x-velocity vert. grd                     |
| duydz_bar         | xc, yc           | 1/s         | Mean y-velocity vert. grd                     |
| taud_acx          | xc, yc           | Pa          | Driving stress (x-dir)                        |
| taud_acy          | xc, yc           | Pa          | Driving stress (y-dir)                        |
| taud              | xc, yc           | Pa          | Driving stress magnitude                      |
| taub_acx          | xc, yc           | Pa          | Basal stress (x-dir)                          |
| taub_acy          | xc, yc           | Pa          | Basal stress (y-dir)                          |
| taub              | xc, yc           | Pa          | Basal stress magnitude                        |
| taul_int_acx      | xc, yc           | Pa          | Internal stress (x-dir)                       |
| taul_int_acy      | xc, yc           | Pa          | Internal stress (y-dir)                       |
| qq_gl_acx         | xc, yc           | m^2/s       | Flux across grounding ln                      |
| qq_gl_acy         | xc, yc           | m^2/s       | Flux across grounding ln                      |
| qq_acx            | xc, yc           | m^2/s       | Flux (x-dir)                                  |
| qq_acy            | xc, yc           | m^2/s       | Flux (y-dir)                                  |
| qq                | xc, yc           | m^2/s       | Flux magnitude                                |
| de_eff            | xc, yc, zeta     | 1/s         | Effective strain rate                         |
| visc_eff          | xc, yc, zeta     | Pa s        | Effective viscosity                           |
| visc_eff_int      | xc, yc           | Pa s        | Depth-integrated viscosity                    |
| N_eff             | xc, yc           | Pa          | Effective pressure                            |
| cb_tgt            | xc, yc           | -           | Target basal parameter                        |
| cb_ref            | xc, yc           | -           | Reference basal parameter                     |
| c_bed             | xc, yc           | -           | Basal drag coefficient                        |
| beta_acx          | xc, yc           | -           | Basal stress factor (x)                       |
| beta_acy          | xc, yc           | -           | Basal stress factor (y)                       |
| beta              | xc, yc           | -           | Basal stress factor mag.                      |
| beta_eff          | xc, yc           | -           | Effective basal factor                        |
| f_vbvs            | xc, yc           | -           | Vertical basal stress                         |
| ssa_mask_acx      | xc, yc           | -           | SSA mask (x-dir)                              |
| ssa_mask_acy      | xc, yc           | -           | SSA mask (y-dir)                              |
| ssa_err_acx       | xc, yc           | -           | SSA error (x-dir)                             |
| ssa_err_acy       | xc, yc           | -           | SSA error (y-dir)                             |
