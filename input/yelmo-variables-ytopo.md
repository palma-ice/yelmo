# ytopo

| variable          | dimensions  | units       | long_name                                          |
|-------------------|-------------|-------------|----------------------------------------------------|
| H_ice             | xc, yc      | m           | Ice thickness                                      |
| dHidt             | xc, yc      | m/yr        | Ice thickness rate of change                       |
| dHidt_dyn         | xc, yc      | m/yr        | Ice thickness change due to dynamics               |
| mb_net            | xc, yc      | m/yr        | Actual mass balance applied                        |
| mb_relax          | xc, yc      | m/yr        | Change in mass balance due to relaxation           |
| mb_resid          | xc, yc      | m/yr        | Residual mass balance                              |
| mb_err            | xc, yc      | m/yr        | Residual error in mass balance accounting          |
| smb               | xc, yc      | m/yr        | Surface mass balance                               |
| bmb               | xc, yc      | m/yr        | Combined basal mass balance                        |
| fmb               | xc, yc      | m/yr        | Combined frontal mass balance                      |
| dmb               | xc, yc      | m/yr        | Subgrid discharge mass balance                     |
| cmb               | xc, yc      | m/yr        | Calving mass balance                               |
| bmb_ref           | xc, yc      | m/yr        | Reference basal mass balance                       |
| fmb_ref           | xc, yc      | m/yr        | Reference frontal mass balance                     |
| dmb_ref           | xc, yc      | m/yr        | Reference subgrid discharge mass balance           |
| cmb_flt           | xc, yc      | m/yr        | Floating calving rate                              |
| cmb_grnd          | xc, yc      | m/yr        | Grounded calving rate                              |
| z_srf             | xc, yc      | m           | Surface elevation                                  |
| dzsdt             | xc, yc      | m/yr        | Surface elevation rate of change                   |
| mask_adv          | xc, yc      |             | Advection mask                                     |
| eps_eff           | xc, yc      | 1/yr        | Effective strain rate                              |
| tau_eff           | xc, yc      | Pa          | Effective stress                                   |
| z_base            | xc, yc      | m           | Ice-base elevation                                 |
| dzsdx             | xc, yc      | m/m         | Surface elevation slope, acx nodes                 |
| dzsdy             | xc, yc      | m/m         | Surface elevation slope, acy nodes                 |
| dHidx             | xc, yc      | m/m         | Ice thickness gradient, acx nodes                  |
| dHidy             | xc, yc      | m/m         | Ice thickness gradient, acy nodes                  |
| dzbdx             | xc, yc      | m/m         | Bedrock slope, acx nodes                           |
| dzbdy             | xc, yc      | m/m         | Bedrock slope, acy nodes                           |
| H_eff             | xc, yc      | m           | Effective ice thickness (margin-corrected)         |
| H_grnd            | xc, yc      | m           | Grounded ice thickness                             |
| H_calv            | xc, yc      | m           | Calving parameter field, ice thickness limit       |
| kt                | xc, yc      |             | Calving parameter field, vm-l19                    |
| z_bed_filt        | xc, yc      | m           | Filtered bedrock elevation                         |
| f_grnd            | xc, yc      |             | Grounded fraction                                  |
| f_grnd_acx        | xc, yc      |             | Grounded fraction (acx nodes)                      |
| f_grnd_acy        | xc, yc      |             | Grounded fraction (acy nodes)                      |
| f_grnd_ab         | xc, yc      |             | Grounded fraction (ab nodes)                       |
| f_ice             | xc, yc      |             | Ice-covered fraction                               |
| f_grnd_bmb        | xc, yc      |             | Grounded fraction for basal mass balance           |
| f_grnd_pin        | xc, yc      |             | Grounded fraction from subgrid pinning points      |
| dist_margin       | xc, yc      | m           | Distance to nearest margin point                   |
| dist_grline       | xc, yc      | m           | Distance to nearest grounding-line point           |
| mask_bed          | xc, yc      |             | Multi-valued bed mask                              |
| mask_grz          | xc, yc      |             | Multi-valued grounding-line zone mask              |
| mask_frnt         | xc, yc      |             | Multi-valued ice front mask                        |
| dHidt_dyn_n       | xc, yc      | m/yr        | Ice thickness change due to advection (previous)   |
| H_ice_n           | xc, yc      | m           | Ice thickness from previous timestep               |
| z_srf_n           | xc, yc      | m           | Surface elevation from previous timestep           |
| H_ice_dyn         | xc, yc      | m           | Dynamic ice thickness                              |
| f_ice_dyn         | xc, yc      |             | Dynamic ice-covered fraction                       |
