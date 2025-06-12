# ytopo

| id | variable          | dimensions  | units       | long_name                                          |
|----|-------------------|-------------|-------------|----------------------------------------------------|
|  1 | H_ice             | xc, yc      | m           | Ice thickness                                      |
|  2 | dHidt             | xc, yc      | m/yr        | Ice thickness rate of change                       |
|  3 | dHidt_dyn         | xc, yc      | m/yr        | Ice thickness change due to dynamics               |
|  4 | mb_net            | xc, yc      | m/yr        | Actual mass balance applied                        |
|  5 | mb_relax          | xc, yc      | m/yr        | Change in mass balance due to relaxation           |
|  6 | mb_resid          | xc, yc      | m/yr        | Residual mass balance                              |
|  7 | mb_err            | xc, yc      | m/yr        | Residual error in mass balance accounting          |
|  8 | smb               | xc, yc      | m/yr        | Surface mass balance                               |
|  9 | bmb               | xc, yc      | m/yr        | Combined basal mass balance                        |
| 10 | fmb               | xc, yc      | m/yr        | Combined frontal mass balance                      |
| 11 | dmb               | xc, yc      | m/yr        | Subgrid discharge mass balance                     |
| 12 | cmb               | xc, yc      | m/yr        | Calving mass balance                               |
| 13 | bmb_ref           | xc, yc      | m/yr        | Reference basal mass balance                       |
| 14 | fmb_ref           | xc, yc      | m/yr        | Reference frontal mass balance                     |
| 15 | dmb_ref           | xc, yc      | m/yr        | Reference subgrid discharge mass balance           |
| 16 | cmb_flt           | xc, yc      | m/yr        | Floating calving rate                              |
| 17 | cmb_grnd          | xc, yc      | m/yr        | Grounded calving rate                              |
| 18 | z_srf             | xc, yc      | m           | Surface elevation                                  |
| 19 | dzsdt             | xc, yc      | m/yr        | Surface elevation rate of change                   |
| 20 | mask_adv          | xc, yc      |             | Advection mask                                     |
| 21 | eps_eff           | xc, yc      | 1/yr        | Effective strain rate                              |
| 22 | tau_eff           | xc, yc      | Pa          | Effective stress                                   |
| 23 | z_base            | xc, yc      | m           | Ice-base elevation                                 |
| 24 | dzsdx             | xc, yc      | m/m         | Surface elevation slope, acx nodes                 |
| 25 | dzsdy             | xc, yc      | m/m         | Surface elevation slope, acy nodes                 |
| 26 | dHidx             | xc, yc      | m/m         | Ice thickness gradient, acx nodes                  |
| 27 | dHidy             | xc, yc      | m/m         | Ice thickness gradient, acy nodes                  |
| 28 | dzbdx             | xc, yc      | m/m         | Bedrock slope, acx nodes                           |
| 29 | dzbdy             | xc, yc      | m/m         | Bedrock slope, acy nodes                           |
| 30 | H_eff             | xc, yc      | m           | Effective ice thickness (margin-corrected)         |
| 31 | H_grnd            | xc, yc      | m           | Grounded ice thickness                             |
| 32 | H_calv            | xc, yc      | m           | Calving parameter field, ice thickness limit       |
| 33 | kt_calv           | xc, yc      |             | Calving parameter field, vm-l19                    |
| 34 | z_bed_filt        | xc, yc      | m           | Filtered bedrock elevation                         |
| 35 | f_grnd            | xc, yc      |             | Grounded fraction                                  |
| 36 | f_grnd_acx        | xc, yc      |             | Grounded fraction (acx nodes)                      |
| 37 | f_grnd_acy        | xc, yc      |             | Grounded fraction (acy nodes)                      |
| 38 | f_grnd_ab         | xc, yc      |             | Grounded fraction (ab nodes)                       |
| 39 | f_ice             | xc, yc      |             | Ice-covered fraction                               |
| 40 | f_grnd_bmb        | xc, yc      |             | Grounded fraction for basal mass balance           |
| 41 | f_grnd_pin        | xc, yc      |             | Grounded fraction from subgrid pinning points      |
| 42 | dist_margin       | xc, yc      | m           | Distance to nearest margin point                   |
| 43 | dist_grline       | xc, yc      | m           | Distance to nearest grounding-line point           |
| 44 | mask_bed          | xc, yc      |             | Multi-valued bed mask                              |
| 45 | mask_grz          | xc, yc      |             | Multi-valued grounding-line zone mask              |
| 46 | mask_frnt         | xc, yc      |             | Multi-valued ice front mask                        |
| 47 | dHidt_dyn_n       | xc, yc      | m/yr        | Ice thickness change due to advection (previous)   |
| 48 | H_ice_n           | xc, yc      | m           | Ice thickness from previous timestep               |
| 49 | z_srf_n           | xc, yc      | m           | Surface elevation from previous timestep           |
| 50 | H_ice_dyn         | xc, yc      | m           | Dynamic ice thickness                              |
| 51 | f_ice_dyn         | xc, yc      |             | Dynamic ice-covered fraction                       |
| 52 | pc_pred_H_ice     | xc, yc      | m           | Predicted ice thickness                            |
| 53 | pc_pred_dHidt_dyn | xc, yc      | m/yr        | Predicted dynamic ice thickness rate of change     |
| 54 | pc_pred_mb_net    | xc, yc      | m/yr        | Predicted net mass balance                         |
| 55 | pc_pred_smb       | xc, yc      | m/yr        | Predicted surface mass balance                     |
| 56 | pc_pred_bmb       | xc, yc      | m/yr        | Predicted basal mass balance                       |
| 57 | pc_pred_fmb       | xc, yc      | m/yr        | Predicted frontal mass balance                     |
| 58 | pc_pred_dmb       | xc, yc      | m/yr        | Predicted discharge mass balance                   |
| 59 | pc_pred_cmb       | xc, yc      | m/yr        | Predicted calving mass balance                     |
| 60 | pc_corr_H_ice     | xc, yc      | m           | Corrected ice thickness                            |
| 61 | pc_corr_dHidt_dyn | xc, yc      | m/yr        | Corrected dynamic ice thickness rate of change     |
| 62 | pc_corr_mb_net    | xc, yc      | m/yr        | Corrected net mass balance                         |
| 63 | pc_corr_smb       | xc, yc      | m/yr        | Corrected surface mass balance                     |
| 64 | pc_corr_bmb       | xc, yc      | m/yr        | Corrected basal mass balance                       |
| 65 | pc_corr_fmb       | xc, yc      | m/yr        | Corrected frontal mass balance                     |
| 66 | pc_corr_dmb       | xc, yc      | m/yr        | Corrected discharge mass balance                   |
| 67 | pc_corr_cmb       | xc, yc      | m/yr        | Corrected calving mass balance                     |
| 68 | lsf               | xc, yc      |             | Level-set function                                 |
| 69 | dlsfdt            | xc, yc      |             | difference Level-set function                      |
| 70 | cmb_flt_x         | xc, yc      | m/yr        | Floating calving rate (x-direction)                |
| 71 | cmb_flt_y         | xc, yc      | m/yr        | Floating calving rate (y-direction)                |
