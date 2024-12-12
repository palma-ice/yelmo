# ymat

| variable          | dimensions       | units       | long_name                                     |
|-------------------|------------------|-------------|-----------------------------------------------|
| enh               | xc, yc, zeta     | -           | Enhancement factor                            |
| enh_bnd           | xc, yc, zeta     | -           | Imposed enhancement factor                    |
| enh_bar           | xc, yc           | -           | Depth-averaged enhancement                    |
| ATT               | xc, yc, zeta     | -           | Rate factor                                   |
| ATT_bar           | xc, yc           | -           | Depth-averaged rate factor                    |
| visc              | xc, yc, zeta     | Pa·yr       | Ice viscosity                                 |
| visc_bar          | xc, yc           | Pa·yr       | Depth-averaged ice viscosity                  |
| visc_int          | xc, yc           | Pa yr m     | Ice viscosity interpolated at interfaces      |
| f_shear_bar       | xc, yc           | -           | Depth-averaged shear fraction                 |
| dep_time          | xc, yc, zeta     | yr          | Ice deposition time (for online age tracing)  |
| depth_iso         | xc, yc, zeta     | m           | Depth of specific isochronal layers           |
| strn2D_dxx        | xc, yc           | 1/yr        | 2D strain rate tensor component dxx           |
| strn2D_dyy        | xc, yc           | 1/yr        | 2D strain rate tensor component dyy           |
| strn2D_dxy        | xc, yc           | 1/yr        | 2D strain rate tensor component dxy           |
| strn2D_dxz        | xc, yc           | 1/yr        | 2D strain rate tensor component dxz           |
| strn2D_dyz        | xc, yc           | 1/yr        | 2D strain rate tensor component dyz           |
| strn2D_de         | xc, yc           | 1/yr        | 2D effective strain rate                      |
| strn2D_div        | xc, yc           | 1/yr        | 2D horizontal divergence                      |
| strn2D_f_shear    | xc, yc           |             | 2D strain rate shear fraction                 |
| strn_dxx          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dxx              |
| strn_dyy          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dyy              |
| strn_dxy          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dxy              |
| strn_dxz          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dxz              |
| strn_dyz          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dyz              |
| strn_de           | xc, yc, zeta     | 1/yr        | Effective strain rate                         |
| strn_div          | xc, yc, zeta     | 1/yr        | Horizontal divergence                         |
| strn_f_shear      | xc, yc, zeta     |             | Strain rate shear fraction                    |
| strs2D_txx        | xc, yc           | Pa          | 2D stress tensor component txx                |
| strs2D_tyy        | xc, yc           | Pa          | 2D stress tensor component tyy                |
| strs2D_txy        | xc, yc           | Pa          | 2D stress tensor component txy                |
| strs2D_txz        | xc, yc           | Pa          | 2D stress tensor component txz                |
| strs2D_tyz        | xc, yc           | Pa          | 2D stress tensor component tyz                |
| strs2D_te         | xc, yc           | Pa          | 2D effective stress                           |
| strs2D_tau_eig_1  | xc, yc           | Pa          | 2D stress first principal eigenvalue          |
| strs2D_tau_eig_2  | xc, yc           | Pa          | 2D stress second principal eigenvalue         |
| strs_txx          | xc, yc, zeta     | Pa          | Stress tensor component txx                   |
| strs_tyy          | xc, yc, zeta     | Pa          | Stress tensor component tyy                   |
| strs_txy          | xc, yc, zeta     | Pa          | Stress tensor component txy                   |
| strs_txz          | xc, yc, zeta     | Pa          | Stress tensor component txz                   |
| strs_tyz          | xc, yc, zeta     | Pa          | Stress tensor component tyz                   |
| strs_te           | xc, yc, zeta     | Pa          | Effective stress                              |
