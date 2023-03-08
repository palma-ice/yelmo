{
    "defaults" :
    {
        "jobname"       : "Yelmo",
        "email"         : "USER@ucm.es",
        "group"         : "anthroia", 
        "omp"           : 0,
        "wall"          : 24, 
        "qos"           : "priority",
        "partition"     : "haswell",
        "job_template"  : "config/pik_submit_slurm"
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
