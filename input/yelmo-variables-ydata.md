# ydata

| variable          | dimensions    | units       | long_name                                          |
|-------------------|---------------|-------------|----------------------------------------------------|
| pd_H_ice          | xc, yc        | m           | Data ice thickness                                 |
| pd_z_srf          | xc, yc        | m           | Data surface elevation                             |
| pd_z_bed          | xc, yc        | m           | Data bedrock elevation                             |
| pd_H_grnd         | xc, yc        | m           | Data overburden ice thickness                      |
| pd_mask_bed       | xc, yc        | -           | Data mask                                          |
| pd_ux_s           | xc, yc        | m/yr        | Data surface velocity in the x-direction           |
| pd_uy_s           | xc, yc        | m/yr        | Data surface velocity in the y-direction           |
| pd_uxy_s          | xc, yc        | m/yr        | Data surface velocity magnitude                    |
| pd_T_srf          | xc, yc        | K           | Data surface temperature                           |
| pd_smb_ref        | xc, yc        | m/yr        | Data surface mass balance                          |
| pd_depth_iso      | xc, yc, zeta  | m           | Data depth of specific isochrones                  |
| err_H_ice         | xc, yc        | m           | Data error in ice thickness                        |
| err_z_srf         | xc, yc        | m           | Data error in surface elevation                    |
| err_z_bed         | xc, yc        | m           | Data error in bedrock elevation                    |
| err_uxy_s         | xc, yc        | m/yr        | Data error in surface velocity magnitude           |
| err_depth_iso     | xc, yc, zeta  | m           | Data error in isochrone depth                      |
