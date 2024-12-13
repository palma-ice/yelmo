# ybound

| id | variable          | dimensions  | units       | long_name                                          |
|----|-------------------|-------------|-------------|----------------------------------------------------|
|  1 | z_bed             | xc, yc      | m           | Bedrock elevation                                  |
|  2 | z_bed_sd          | xc, yc      | m           | Standard deviation of bedrock elevation            |
|  3 | z_sl              | xc, yc      | m           | Sea level elevation                                |
|  4 | H_sed             | xc, yc      | m           | Sediment thickness                                 |
|  5 | smb_ref           | xc, yc      | m/yr        | Surface mass balance                               |
|  6 | T_srf             | xc, yc      | K           | Surface temperature                                |
|  7 | bmb_shlf          | xc, yc      | m/yr        | Basal mass balance for ice shelf                   |
|  8 | fmb_shlf          | xc, yc      | m/yr        | Frontal mass balance for ice shelf                 |
|  9 | T_shlf            | xc, yc      | K           | Ice shelf temperature                              |
| 10 | Q_geo             | xc, yc      | mW m^-2     | Geothermal heat flow at depth                      |
| 11 | enh_srf           | xc, yc      | -           | Enhancement factor at the surface                  |
| 12 | basins            | xc, yc      | -           | Basin identification numbers                       |
| 13 | basin_mask        | xc, yc      | -           | Mask for basins                                    |
| 14 | regions           | xc, yc      | -           | Region identification numbers                      |
| 15 | region_mask       | xc, yc      | -           | Mask for regions                                   |
| 16 | ice_allowed       | xc, yc      | -           | Locations where ice thickness can be greater than zero |
| 17 | calv_mask         | xc, yc      | -           | Locations where calving is not allowed             |
| 18 | H_ice_ref         | xc, yc      | m           | Reference ice thickness for relaxation routines    |
| 19 | z_bed_ref         | xc, yc      | m           | Reference bedrock elevation for relaxation routines |
