{
    "defaults" :
    {
        "jobname"       : "Yelmo",
        "email"         : "USER@ucm.es",
        "group"         : "", 
        "omp"           : 0,
        "wall"          : 24, 
        "qos"           : "normal",
        "partition"     : "",
        "job_template"  : "config/brigit_submit_slurm"
    },

    "exe_aliases" : 
        {   "benchmarks" : "libyelmo/bin/yelmo_benchmarks.x",
            "mismip"     : "libyelmo/bin/yelmo_mismip.x",
            "initmip"    : "libyelmo/bin/yelmo_initmip.x",
            "opt"        : "libyelmo/bin/yelmo_opt.x",
            "trough"     : "libyelmo/bin/yelmo_trough.x",
            "ismiphom"   : "libyelmo/bin/yelmo_ismiphom.x",
            "slab"       : "libyelmo/bin/yelmo_slab.x",
            "regridding" : "libyelmo/bin/yelmo_test_regridding.x"
        },

    "grp_aliases" : {},

    "par_paths" : {},

    "files" : [], 

    "dir-special" : {},

    "links" : 
        ["input","ice_data"],

    "const_paths" : 
        {   "EISMINT"  : "par/yelmo_const_EISMINT.nml",
            "HALFAR"   : "par/yelmo_const_EISMINT.nml",
            "dome"     : "par/yelmo_const_EISMINT.nml",
            "MISMIP3D" : "par/yelmo_const_MISMIP3D.nml",
            "TROUGH"   : "par/yelmo_const_TROUGH.nml",
            "MISMIP+"  : "par/yelmo_const_TROUGH.nml",
            "ISMIPHOM" : "par/yelmo_const_EISMINT.nml", 
            "slab"     : "par/yelmo_const_EISMINT.nml" 
        },

    "const_path_default" : "par/yelmo_const_Earth.nml",
    
    "job_queues" :
        {  "normal" :
            {   "wall" : 10000 }
        }
}
