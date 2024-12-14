# ymat

| id | variable          | dimensions       | units       | long_name                                     |
|----|-------------------|------------------|-------------|-----------------------------------------------|
|  1 | enh               | xc, yc, zeta     | -           | Enhancement factor                            |
|  2 | enh_bnd           | xc, yc, zeta     | -           | Imposed enhancement factor                    |
|  3 | enh_bar           | xc, yc           | -           | Depth-averaged enhancement                    |
|  4 | ATT               | xc, yc, zeta     | -           | Rate factor                                   |
|  5 | ATT_bar           | xc, yc           | -           | Depth-averaged rate factor                    |
|  6 | visc              | xc, yc, zeta     | Pa yr       | Ice viscosity                                 |
|  7 | visc_bar          | xc, yc           | Pa yr       | Depth-averaged ice viscosity                  |
|  8 | visc_int          | xc, yc           | Pa yr m     | Ice viscosity interpolated at interfaces      |
|  9 | f_shear_bar       | xc, yc           | -           | Depth-averaged shear fraction                 |
| 10 | dep_time          | xc, yc, zeta     | yr          | Ice deposition time (for online age tracing)  |
| 11 | depth_iso         | xc, yc, age_iso  | m           | Depth of specific isochronal layers           |
| 12 | strn2D_dxx        | xc, yc           | 1/yr        | 2D strain rate tensor component dxx           |
| 13 | strn2D_dyy        | xc, yc           | 1/yr        | 2D strain rate tensor component dyy           |
| 14 | strn2D_dxy        | xc, yc           | 1/yr        | 2D strain rate tensor component dxy           |
| 15 | strn2D_dxz        | xc, yc           | 1/yr        | 2D strain rate tensor component dxz           |
| 16 | strn2D_dyz        | xc, yc           | 1/yr        | 2D strain rate tensor component dyz           |
| 17 | strn2D_de         | xc, yc           | 1/yr        | 2D effective strain rate                      |
| 18 | strn2D_div        | xc, yc           | 1/yr        | 2D horizontal divergence                      |
| 19 | strn2D_f_shear    | xc, yc           |             | 2D strain rate shear fraction                 |
| 20 | strn_dxx          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dxx              |
| 21 | strn_dyy          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dyy              |
| 22 | strn_dxy          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dxy              |
| 23 | strn_dxz          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dxz              |
| 24 | strn_dyz          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dyz              |
| 25 | strn_de           | xc, yc, zeta     | 1/yr        | Effective strain rate                         |
| 26 | strn_div          | xc, yc, zeta     | 1/yr        | Horizontal divergence                         |
| 27 | strn_f_shear      | xc, yc, zeta     |             | Strain rate shear fraction                    |
| 28 | strs2D_txx        | xc, yc           | Pa          | 2D stress tensor component txx                |
| 29 | strs2D_tyy        | xc, yc           | Pa          | 2D stress tensor component tyy                |
| 30 | strs2D_txy        | xc, yc           | Pa          | 2D stress tensor component txy                |
| 31 | strs2D_txz        | xc, yc           | Pa          | 2D stress tensor component txz                |
| 32 | strs2D_tyz        | xc, yc           | Pa          | 2D stress tensor component tyz                |
| 33 | strs2D_te         | xc, yc           | Pa          | 2D effective stress                           |
| 34 | strs2D_tau_eig_1  | xc, yc           | Pa          | 2D stress first principal eigenvalue          |
| 35 | strs2D_tau_eig_2  | xc, yc           | Pa          | 2D stress second principal eigenvalue         |
| 36 | strs_txx          | xc, yc, zeta     | Pa          | Stress tensor component txx                   |
| 37 | strs_tyy          | xc, yc, zeta     | Pa          | Stress tensor component tyy                   |
| 38 | strs_txy          | xc, yc, zeta     | Pa          | Stress tensor component txy                   |
| 39 | strs_txz          | xc, yc, zeta     | Pa          | Stress tensor component txz                   |
| 40 | strs_tyz          | xc, yc, zeta     | Pa          | Stress tensor component tyz                   |
| 41 | strs_te           | xc, yc, zeta     | Pa          | Effective stress                              |
