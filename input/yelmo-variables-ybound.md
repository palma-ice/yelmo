# ybound

| variable          | dimensions  | units       | long_name                                          |
|-------------------|-------------|-------------|----------------------------------------------------|
| z_bed             | 2D          | m           | Bedrock elevation                                  |
| z_bed_sd          | 2D          | m           | Standard deviation of bedrock elevation            |
| z_sl              | 2D          | m           | Sea level elevation                                |
| H_sed             | 2D          | m           | Sediment thickness                                 |
| smb               | 2D          | m/yr        | Surface mass balance                               |
| T_srf             | 2D          | K           | Surface temperature                                |
| bmb_shlf          | 2D          | m/yr        | Basal mass balance for ice shelf                   |
| fmb_shlf          | 2D          | m/yr        | Frontal mass balance for ice shelf                 |
| T_shlf            | 2D          | K           | Ice shelf temperature                              |
| Q_geo             | 2D          | W m^-2      | Geothermal heat flow at depth                      |
| enh_srf           | 2D          | -           | Enhancement factor at the surface                  |
| basins            | 2D          | -           | Basin identification numbers                       |
| basin_mask        | 2D          | -           | Mask for basins                                    |
| regions           | 2D          | -           | Region identification numbers                      |
| region_mask       | 2D          | -           | Mask for regions                                   |
| ice_allowed       | 2D          | -           | Locations where ice thickness can be greater than zero |
| calv_mask         | 2D          | -           | Locations where calving is not allowed             |
| H_ice_ref         | 2D          | m           | Reference ice thickness for relaxation routines    |
| z_bed_ref         | 2D          | m           | Reference bedrock elevation for relaxation routines |
