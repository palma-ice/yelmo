{
    "defaults" :
        {
        "jobname"       : "Yelmo",
        "email"         : "USER@awi.de",
        "group"         : "envi",
        "omp"           : 0,
        "wall"          : "48:00:00",
        "qos"           : "48h",
        "partition"     : "haswell",
        "job_template"  : "config/albedo_submit_slurm"
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
        ["input","ice_data","maps"],

    "job_queues" :
        {   "priority" :
            {   "wall" : 24  },
            "short" :
            {   "wall" : 24  },
            "medium" :
            {   "wall" : 168 },
            "long" :
            {   "wall" : 720 }
        }
}
