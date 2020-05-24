fldr=tmp/iso
job run --shell -f -o ${fldr}/test1 -p marine_shelf.obs_f=10.0 snapclim_hybrid.f_to=0.4 snapclim_hybrid.f_hol=0.6,0.8,1.0 snapclim_hybrid.f_glac=1.0 smbpal_par.abl_method=\'pdd\' smbpal_par.mm_snow=2.0,3.0,4.0 -- python run_yelmo.py -x -s -e iso {} par/yelmo_Greenland_iso.nml

# ens1
job sample snapclim.f_p=U?0.04,0.07 snapclim_hybrid.f_glac=U?0.7,1.1 snapclim_hybrid.f_hol=U?0.7,1.0 smbpal_par.mm_snow=U?2.0,4.0 smbpal_par.mm_ice=U?6.0,8.0 --size 100 --seed 4 > ens1.txt
job run --shell -f -o ${fldr}/test11 -i ens1.txt -- python run_yelmo.py -x -s -e iso {} par/yelmo_Greenland_iso.nml

# ens2
job sample snapclim.f_p=U?0.04,0.07 snapclim_hybrid.f_hol=U?0.7,1.0 smbpal_par.mm_snow=U?1.0,4.0 smbpal_par.mm_ice=U?5.0,8.0 itm_par.Pmaxfrac=U?0.5,0.7 --size 100 --seed 4 > ens2.txt
job run --shell -f -o ${fldr}/ens2 -i ens2.txt -- python run_yelmo.py -x -s -e iso {} par/yelmo_Greenland_iso.nml

# ens8
job sample snapclim.f_p=U?0.04,0.07 snapclim_hybrid.f_hol=U?0.7,1.0 smbpal_par.mm_snow=U?1.0,4.0 smbpal_par.mm_ice=U?5.0,8.0 ctrl.dpr_hol=U?0.0,0.5 --size 100 --seed 4 > ens8.txt
job run --shell -f -o ${fldr}/ens8 -i ens8.txt -- python run_yelmo.py -x -s -e iso {} par/yelmo_Greenland_iso.nml
mv ens8.txt ${fldr}/ens8/ens_params.txt 

# Check results
for D in tmp/iso/ens2/* ; do ./check_sim.x $D ; done > tmp/iso/ens2/ens_err.txt