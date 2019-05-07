# yelmo
Yelmo ice sheet model code base



# Getting started 




# Installing NetCDF



# Installing LIS

1. Download the LIS source:
https://www.ssisc.org/lis/

2. Configure the package and install it
in the location of your choice (make sure to enable the Fortran90 interface):
```
cd lis-2.0.18
./configure --prefix=$HOME/apps/lis/2.0.18 --enable-f90
make
make install
make install check

# Make a handy link to this version of LIS
ln -s 2.0.18 lis
```

3. Add LIS path to the LD\_LIBRARY\_PATH in .bash\_profile, .bashrc or .bash\_aliases:
```
# lis library paths
LD_LIBRARY_PATH=$HOME/apps/lis/lis/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
```

That' it. LIS should now be available to use with Yelmo.


