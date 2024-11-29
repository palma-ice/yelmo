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
| dep_time          | xc, yc, zeta     | years       | Ice deposition time (for online age tracing)  |
| depth_iso         | xc, yc, zeta     | m           | Depth of specific isochronal layers           |
