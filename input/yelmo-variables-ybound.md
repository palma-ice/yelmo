# ybound

| variable          | dimensions  | units       | long_name                                          |
|-------------------|-------------|-------------|----------------------------------------------------|
| z_bed             | xc, yc      | m           | Bedrock elevation                                  |
| z_bed_sd          | xc, yc      | m           | Standard deviation of bedrock elevation            |
| z_sl              | xc, yc      | m           | Sea level elevation                                |
| H_sed             | xc, yc      | m           | Sediment thickness                                 |
| smb               | xc, yc      | m/yr        | Surface mass balance                               |
| T_srf             | xc, yc      | K           | Surface temperature                                |
| bmb_shlf          | xc, yc      | m/yr        | Basal mass balance for ice shelf                   |
| fmb_shlf          | xc, yc      | m/yr        | Frontal mass balance for ice shelf                 |
| T_shlf            | xc, yc      | K           | Ice shelf temperature                              |
| Q_geo             | xc, yc      | W m^-2      | Geothermal heat flow at depth                      |
| enh_srf           | xc, yc      | -           | Enhancement factor at the surface                  |
| basins            | xc, yc      | -           | Basin identification numbers                       |
| basin_mask        | xc, yc      | -           | Mask for basins                                    |
| regions           | xc, yc      | -           | Region identification numbers                      |
| region_mask       | xc, yc      | -           | Mask for regions                                   |
| ice_allowed       | xc, yc      | -           | Locations where ice thickness can be greater than zero |
| calv_mask         | xc, yc      | -           | Locations where calving is not allowed             |
| H_ice_ref         | xc, yc      | m           | Reference ice thickness for relaxation routines    |
| z_bed_ref         | xc, yc      | m           | Reference bedrock elevation for relaxation routines |
