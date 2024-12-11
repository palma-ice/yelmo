# ydata

| variable          | dimensions    | units       | long_name                                          |
|-------------------|---------------|-------------|----------------------------------------------------|
| H_ice             | xc, yc        | m           | Data ice thickness                                 |
| z_srf             | xc, yc        | m           | Data surface elevation                             |
| z_bed             | xc, yc        | m           | Data bedrock elevation                             |
| H_grnd            | xc, yc        | m           | Data overburden ice thickness                      |
| mask_bed          | xc, yc        | -           | Data mask                                          |
| ux_s              | xc, yc        | m/yr        | Data surface velocity in the x-direction           |
| uy_s              | xc, yc        | m/yr        | Data surface velocity in the y-direction           |
| uxy_s             | xc, yc        | m/yr        | Data surface velocity magnitude                    |
| T_srf             | xc, yc        | K           | Data surface temperature                           |
| smb               | xc, yc        | m/yr        | Data surface mass balance                          |
| depth_iso         | xc, yc, zeta  | m           | Data depth of specific isochrones                  |
| err_H_ice         | xc, yc        | m           | Data error in ice thickness                        |
| err_z_srf         | xc, yc        | m           | Data error in surface elevation                    |
| err_z_bed         | xc, yc        | m           | Data error in bedrock elevation                    |
| err_uxy_s         | xc, yc        | m/yr        | Data error in surface velocity magnitude           |
| err_depth_iso     | xc, yc, zeta  | m           | Data error in isochrone depth                      |
