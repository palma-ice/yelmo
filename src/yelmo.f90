module yelmo
    ! Wrapper to hold all modules needed for libyelmo 

    use yelmo_defs 
    use yelmo_grid, only : yelmo_init_grid, yelmo_grid_write
    use yelmo_io 
    use yelmo_ice 

!     use yelmo_timesteps 
    
!     use yelmo_topography
!     use yelmo_boundaries 
    use yelmo_regions, only : yelmo_write_reg_init, yelmo_write_reg_step 
    
    use basal_dragging, only : calc_cb_ref 
    
    use topography, only : mask_bed_ocean, mask_bed_land, mask_bed_frozen, &
                    mask_bed_stream, mask_bed_grline, mask_bed_float , mask_bed_island
    
end module yelmo
