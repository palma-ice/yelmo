module yelmo
    ! Wrapper to hold all modules needed for libyelmo 

    use yelmo_defs 
    use yelmo_grid, only : yelmo_init_grid, yelmo_grid_write
    use yelmo_io 
    use yelmo_ice 

!     use yelmo_timesteps 
    
!     use yelmo_topography
!     use yelmo_boundaries 
    use yelmo_regions, only : write_yreg_init, write_yreg_step 
    
end module yelmo 
