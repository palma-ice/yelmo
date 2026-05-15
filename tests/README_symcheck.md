# SSA-energy reflection-symmetry checks

Two layers of testing protect against the L↔R / T↔B asymmetry bug that was
fixed in `solver_ssa_ac_energy`: a fast unit test on the linear assembler
itself, and a slower run-based check on the integrated H_ice field of a
known-symmetric setup (CalvingMIP exp1 on its circular domain).

## What this is checking

The energy SSA assembler builds `K · u = b` over the C-grid. For an input
configuration that is mirror-symmetric about a vertical (or horizontal) line,
both `K` and `b` must be equivariant under the C-grid reflection `P`:

```
ux(i, j) at right face of cell i  ->  ux(nx-i,   j) with sign -1
uy(i, j) at top face of cell i    ->  uy(nx+1-i, j) with sign +1
```

If `K` and `b` respect this symmetry, the linear solve produces a
reflection-symmetric `u`, and any downstream physics (viscosity update,
mass transport, …) that is itself symmetry-preserving will keep H_ice
mirror-symmetric.

`test_ssa_energy_sym` (existing) checks `K = Kᵀ`, which CG needs but is
**not** the same property. `test_ssa_energy_lr_sym` (here) is the missing
reflection-equivariance check that catches the bug `test_ssa_energy_sym`
let through.

## Layer 1 — unit test on the assembler

Direct, deterministic, no I/O. Builds three mirror-symmetric ice
configurations (vertical strip, horizontal strip, square patch), turns on
`taul_int` at the calving fronts, calls
`linear_solver_matrix_ssa_ac_csr_2D_energy` once, and walks the assembled
`K` and `b` for reflection equivariance. Fails fast on regression.

```bash
# Build and run (from repo root):
make ssa_energy_lr_sym
./libyelmo/bin/test_ssa_energy_lr_sym.x
echo $?
```

Expected output ends with:

```
 ALL CASES PASS: K and b are reflection-symmetric
```

and exit status 0. On regression the program prints the failing rows
(`b x row=…`, `K x (n,c)=…`) with the offending values, and exits with
status 1.

## Layer 2 — run-based check on H_ice (CalvingMIP exp1)

The unit test only exercises the assembler. To catch regressions further
downstream (e.g. a non-equivariant viscosity stencil, an asymmetric LSF
advection step, a beta-staggering bug at the grounding line), run a short
CalvingMIP exp1 with both `ssa_solver = "energy"` and `ssa_solver =
"residual"` and compare H_ice symmetry:

```bash
# 1. Stage two run dirs differing only in ssa_solver.
WT=$(pwd)
for s in energy residual; do
    mkdir -p tmp/sym/$s/exp1
    cp par/yelmo_calvingmip.nml tmp/sym/$s/exp1/
    ln -sf $WT/input    tmp/sym/$s/exp1/input
    ln -sf $WT/maps     tmp/sym/$s/exp1/maps
    ln -sf $WT/ice_data tmp/sym/$s/exp1/ice_data
done
sed -i'.bak' 's/ssa_solver *= *"energy"/ssa_solver = "residual"/' tmp/sym/residual/exp1/yelmo_calvingmip.nml

# 2. Run both. CalvingMIP exp1 default is 10000 yr; for a quick sanity
#    check shorten time_end by editing the nml. ~17 min energy, ~10 s
#    residual on a current Mac.
make calving
for s in energy residual; do
    (cd tmp/sym/$s/exp1 && $WT/libyelmo/bin/yelmo_calving.x yelmo_calvingmip.nml > run.log 2>&1)
done

# 3. Diff symmetry. The residual baseline reaches machine precision
#    (~1e-7); the energy solver should match it to within a few percent.
julia tests/symcheck.jl tmp/sym/energy/exp1/yelmo2D.nc tmp/sym/residual/exp1/yelmo2D.nc
```

Healthy output (post-fix):

```
tmp/sym/energy/exp1/yelmo2D.nc
  shape=(65, 65)  Hmax= 1296.22  …
  asym L-R   : Linf=4.492e-03  L1=2.553e-04
  …
tmp/sym/residual/exp1/yelmo2D.nc
  shape=(65, 65)  Hmax= 1281.16  …
  asym L-R   : Linf=2.382e-07  L1=3.154e-08
  …
```

Regression output (pre-fix, for reference):

```
  asym L-R   : Linf=7.189e-01  L1=1.817e-01    # 72% of Hmax
```

`Hmax` for the energy solver should land within ~1% of the residual baseline
(~1281 m for the 25 km / 10 kyr exp1 setup); a much larger value indicates
the front stress is being mis-applied.

### Julia dependencies

`tests/symcheck.jl` needs `NCDatasets` (everything else is stdlib). One-time
setup in your global Julia env:

```julia
import Pkg; Pkg.add("NCDatasets")
```

If you prefer Python, the original `symcheck.py` predecessor lived in
`tmp/symcheck.py` during development; it does the same thing with `netCDF4`.
