module yelmo_c_api
  
  use yelmo

  implicit none

  type(yelmo_class), target :: ylmo1   ! persistent state
  type(yelmo_class), target :: ylmo2   ! persistent state
  type(yelmo_class), pointer :: ylmo    ! persistent state
  
  ! type(yelmo_class) :: ylmo    ! persistent state

contains

  subroutine yelmo_set_alias(alias) bind(C, name="yelmo_set_alias")
    use iso_c_binding
    character(c_char), intent(in) :: alias(*)

    character(len=56)             :: f_alias

    f_alias = trim(c_to_f_string(alias))

    if (trim(f_alias) .eq. "ylmo2") then
      ylmo => ylmo2
    else
      ylmo => ylmo1
    end if

    return
  end subroutine yelmo_set_alias

  subroutine yelmo_init_wrapper(filename, grid_def, time, alias) bind(C, name="yelmo_init")
    use iso_c_binding
    character(c_char), intent(in) :: filename(*)
    character(c_char), intent(in) :: grid_def(*)
    real(c_double), value         :: time
    character(c_char), intent(in) :: alias(*)

    ! Local variables
    character(len=1028)           :: f_filename
    character(len=56)             :: f_grid_def

    ! Convert C string to Fortran string (null-terminated scan)
    f_filename    = trim(c_to_f_string(filename))
    f_grid_def    = trim(c_to_f_string(grid_def))
    
    call yelmo_set_alias(alias)

    call yelmo_init(ylmo, filename=trim(f_filename), grid_def=trim(f_grid_def), time=real(time, wp))

    nullify(ylmo)

  end subroutine


  subroutine yelmo_init_state_wrapper(time, thrm_method, alias) bind(C, name="yelmo_init_state")
    use iso_c_binding
    real(c_double), value         :: time
    character(c_char), intent(in) :: thrm_method(*) ! e.g., "robin-cold"
    character(c_char), intent(in) :: alias(*)

    ! Local variables
    character(len=56) :: f_thrm_method

    ! Convert C string to Fortran string (null-terminated scan)
    f_thrm_method = trim(c_to_f_string(thrm_method))

    call yelmo_set_alias(alias)

    call yelmo_init_state(ylmo, time=real(time, wp), thrm_method=trim(f_thrm_method))

    nullify(ylmo)

  end subroutine

  subroutine yelmo_step_wrapper(time,alias) bind(C, name="yelmo_step")
    use iso_c_binding
    real(c_double), value :: time
    character(c_char), intent(in) :: alias(*)

    call yelmo_set_alias(alias)

    call yelmo_update(ylmo, real(time, wp))

    nullify(ylmo)

  end subroutine

  ! =========================================================================
  ! C-bound grid sizes getter: returns nx, ny, nz_aa, nz_ac only.
  ! Call this first to learn buffer sizes before calling yelmo_get_grid_info.
  ! =========================================================================
  subroutine yelmo_get_grid_sizes(nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac, alias) &
                                  bind(C, name="yelmo_get_grid_sizes")
    use iso_c_binding
    integer(c_int), intent(out) :: nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac
    character(c_char), intent(in) :: alias(*)
    
    call yelmo_set_alias(alias)

    nx    = ylmo%grd%nx
    ny    = ylmo%grd%ny
    nz_aa = ylmo%par%nz_aa
    nz_ac = ylmo%par%nz_ac
    nzr_aa = ylmo%thrm%par%nzr_aa
    nzr_ac = ylmo%thrm%par%nzr_ac

    nullify(ylmo)

  end subroutine yelmo_get_grid_sizes


  ! =========================================================================
  ! C-bound grid axes getter: fills pre-allocated buffers with axis values.
  ! Caller must have already called yelmo_get_grid_sizes and allocated:
  !   xc[nx], yc[ny], zeta_aa[nz_aa], zeta_ac[nz_ac]
  ! =========================================================================
  subroutine yelmo_get_grid_info(xc, yc, zeta_aa, zeta_ac, zeta_r_aa, zeta_r_ac, alias) &
                                 bind(C, name="yelmo_get_grid_info")
    use iso_c_binding
    type(c_ptr), value :: xc, yc, zeta_aa, zeta_ac, zeta_r_aa, zeta_r_ac
    character(c_char), intent(in) :: alias(*)

    ! Local variables
    real(c_double), pointer :: xc_f(:)      => null()
    real(c_double), pointer :: yc_f(:)      => null()
    real(c_double), pointer :: zeta_aa_f(:) => null()
    real(c_double), pointer :: zeta_ac_f(:) => null()
    real(c_double), pointer :: zeta_r_aa_f(:) => null()
    real(c_double), pointer :: zeta_r_ac_f(:) => null()

    call yelmo_set_alias(alias)

    if (c_associated(xc)) then
      call c_f_pointer(xc,      xc_f,      [ylmo%grd%nx])
      xc_f = real(ylmo%grd%xc, c_double)
    end if

    if (c_associated(yc)) then
      call c_f_pointer(yc,      yc_f,      [ylmo%grd%ny])
      yc_f = real(ylmo%grd%yc, c_double)
    end if

    if (c_associated(zeta_aa)) then
      call c_f_pointer(zeta_aa, zeta_aa_f, [ylmo%par%nz_aa])
      zeta_aa_f = real(ylmo%par%zeta_aa, c_double)
    end if

    if (c_associated(zeta_ac)) then
      call c_f_pointer(zeta_ac, zeta_ac_f, [ylmo%par%nz_ac])
      zeta_ac_f = real(ylmo%par%zeta_ac, c_double)
    end if

    if (c_associated(zeta_r_aa)) then
      call c_f_pointer(zeta_r_aa, zeta_r_aa_f, [ylmo%thrm%par%nzr_aa])
      zeta_r_aa_f = real(ylmo%thrm%par%zr%zeta_aa, c_double)
    end if

    if (c_associated(zeta_r_ac)) then
      call c_f_pointer(zeta_r_ac, zeta_r_ac_f, [ylmo%thrm%par%nzr_ac])
      zeta_r_ac_f = real(ylmo%thrm%par%zr%zeta_ac, c_double)
    end if

    nullify(ylmo)

  end subroutine yelmo_get_grid_info

  ! =========================================================================
  ! 2D variable getter
  ! =========================================================================
  subroutine yelmo_get_var2D(v2D, nx, ny, name, alias) bind(C, name="yelmo_get_var2D")
    use iso_c_binding
    integer(c_int), value       :: nx, ny
    character(c_char), intent(in) :: name(*)
    real(c_double), intent(out)   :: v2D(nx, ny)
    character(c_char), intent(in) :: alias(*)
    
    ! Local variables
    character(len=56) :: f_name

    ! Convert C string to Fortran string (null-terminated scan)
    f_name = trim(c_to_f_string(name))

    call yelmo_set_alias(alias)

    select case(trim(f_name))

      ! -----------------------------------------------------------------------
      ! ybound
      ! -----------------------------------------------------------------------
      case("bnd_z_bed");         v2D = real(ylmo%bnd%z_bed,         c_double)
      case("bnd_z_bed_sd");      v2D = real(ylmo%bnd%z_bed_sd,      c_double)
      case("bnd_z_sl");          v2D = real(ylmo%bnd%z_sl,          c_double)
      case("bnd_H_sed");         v2D = real(ylmo%bnd%H_sed,         c_double)
      case("bnd_smb");           v2D = real(ylmo%bnd%smb,           c_double)
      case("bnd_T_srf");         v2D = real(ylmo%bnd%T_srf,         c_double)
      case("bnd_bmb_shlf");      v2D = real(ylmo%bnd%bmb_shlf,      c_double)
      case("bnd_fmb_shlf");      v2D = real(ylmo%bnd%fmb_shlf,      c_double)
      case("bnd_T_shlf");        v2D = real(ylmo%bnd%T_shlf,        c_double)
      case("bnd_Q_geo");         v2D = real(ylmo%bnd%Q_geo,         c_double)
      case("bnd_enh_srf");       v2D = real(ylmo%bnd%enh_srf,       c_double)
      case("bnd_basins");        v2D = real(ylmo%bnd%basins,        c_double)
      case("bnd_basin_mask");    v2D = real(ylmo%bnd%basin_mask,    c_double)
      case("bnd_regions");       v2D = real(ylmo%bnd%regions,       c_double)
      case("bnd_region_mask");   v2D = real(ylmo%bnd%region_mask,   c_double)
      case("bnd_H_ice_ref");     v2D = real(ylmo%bnd%H_ice_ref,     c_double)
      case("bnd_z_bed_ref");     v2D = real(ylmo%bnd%z_bed_ref,     c_double)
      case("bnd_tau_relax");     v2D = real(ylmo%bnd%tau_relax,     c_double)
      case("bnd_z_bed_corr");    v2D = real(ylmo%bnd%z_bed_corr,    c_double)
      case("bnd_dzbdt_corr");    v2D = real(ylmo%bnd%dzbdt_corr,    c_double)

      ! -----------------------------------------------------------------------
      ! ytopo%now (ytopo_state_class)
      ! -----------------------------------------------------------------------
      case("tpo_H_ice");         v2D = real(ylmo%tpo%now%H_ice,         c_double)
      case("tpo_dHidt");         v2D = real(ylmo%tpo%now%dHidt,         c_double)
      case("tpo_dHidt_dyn");     v2D = real(ylmo%tpo%now%dHidt_dyn,     c_double)
      case("tpo_mb_net");        v2D = real(ylmo%tpo%now%mb_net,        c_double)
      case("tpo_mb_relax");      v2D = real(ylmo%tpo%now%mb_relax,      c_double)
      case("tpo_mb_resid");      v2D = real(ylmo%tpo%now%mb_resid,      c_double)
      case("tpo_mb_err");        v2D = real(ylmo%tpo%now%mb_err,        c_double)
      case("tpo_smb");           v2D = real(ylmo%tpo%now%smb,           c_double)
      case("tpo_bmb");           v2D = real(ylmo%tpo%now%bmb,           c_double)
      case("tpo_fmb");           v2D = real(ylmo%tpo%now%fmb,           c_double)
      case("tpo_dmb");           v2D = real(ylmo%tpo%now%dmb,           c_double)
      case("tpo_cmb");           v2D = real(ylmo%tpo%now%cmb,           c_double)
      case("tpo_bmb_ref");       v2D = real(ylmo%tpo%now%bmb_ref,       c_double)
      case("tpo_fmb_ref");       v2D = real(ylmo%tpo%now%fmb_ref,       c_double)
      case("tpo_dmb_ref");       v2D = real(ylmo%tpo%now%dmb_ref,       c_double)
      case("tpo_cmb_flt");       v2D = real(ylmo%tpo%now%cmb_flt,       c_double)
      case("tpo_cmb_flt_x");     v2D = real(ylmo%tpo%now%cmb_flt_x,     c_double)
      case("tpo_cmb_flt_y");     v2D = real(ylmo%tpo%now%cmb_flt_y,     c_double)
      case("tpo_cmb_grnd");      v2D = real(ylmo%tpo%now%cmb_grnd,      c_double)
      case("tpo_cmb_grnd_x");    v2D = real(ylmo%tpo%now%cmb_grnd_x,    c_double)
      case("tpo_cmb_grnd_y");    v2D = real(ylmo%tpo%now%cmb_grnd_y,    c_double)
      case("tpo_cr_acx");        v2D = real(ylmo%tpo%now%cr_acx,        c_double)
      case("tpo_cr_acy");        v2D = real(ylmo%tpo%now%cr_acy,        c_double)
      case("tpo_lsf");           v2D = real(ylmo%tpo%now%lsf,           c_double)
      case("tpo_dlsfdt");        v2D = real(ylmo%tpo%now%dlsfdt,        c_double)
      case("tpo_z_srf");         v2D = real(ylmo%tpo%now%z_srf,         c_double)
      case("tpo_dzsdt");         v2D = real(ylmo%tpo%now%dzsdt,         c_double)
      case("tpo_eps_eff");       v2D = real(ylmo%tpo%now%eps_eff,       c_double)
      case("tpo_tau_eff");       v2D = real(ylmo%tpo%now%tau_eff,       c_double)
      case("tpo_z_base");        v2D = real(ylmo%tpo%now%z_base,        c_double)
      case("tpo_dzsdx");         v2D = real(ylmo%tpo%now%dzsdx,         c_double)
      case("tpo_dzsdy");         v2D = real(ylmo%tpo%now%dzsdy,         c_double)
      case("tpo_dHidx");         v2D = real(ylmo%tpo%now%dHidx,         c_double)
      case("tpo_dHidy");         v2D = real(ylmo%tpo%now%dHidy,         c_double)
      case("tpo_dzbdx");         v2D = real(ylmo%tpo%now%dzbdx,         c_double)
      case("tpo_dzbdy");         v2D = real(ylmo%tpo%now%dzbdy,         c_double)
      case("tpo_dzsdx_aa");      v2D = real(ylmo%tpo%now%dzsdx_aa,      c_double)
      case("tpo_dzsdy_aa");      v2D = real(ylmo%tpo%now%dzsdy_aa,      c_double)
      case("tpo_dHidx_aa");      v2D = real(ylmo%tpo%now%dHidx_aa,      c_double)
      case("tpo_dHidy_aa");      v2D = real(ylmo%tpo%now%dHidy_aa,      c_double)
      case("tpo_dzbdx_aa");      v2D = real(ylmo%tpo%now%dzbdx_aa,      c_double)
      case("tpo_dzbdy_aa");      v2D = real(ylmo%tpo%now%dzbdy_aa,      c_double)
      case("tpo_H_eff");         v2D = real(ylmo%tpo%now%H_eff,         c_double)
      case("tpo_H_grnd");        v2D = real(ylmo%tpo%now%H_grnd,        c_double)
      case("tpo_H_calv");        v2D = real(ylmo%tpo%now%H_calv,        c_double)
      case("tpo_kt");            v2D = real(ylmo%tpo%now%kt,            c_double)
      case("tpo_z_bed_filt");    v2D = real(ylmo%tpo%now%z_bed_filt,    c_double)
      case("tpo_f_grnd");        v2D = real(ylmo%tpo%now%f_grnd,        c_double)
      case("tpo_f_grnd_acx");    v2D = real(ylmo%tpo%now%f_grnd_acx,    c_double)
      case("tpo_f_grnd_acy");    v2D = real(ylmo%tpo%now%f_grnd_acy,    c_double)
      case("tpo_f_grnd_ab");     v2D = real(ylmo%tpo%now%f_grnd_ab,     c_double)
      case("tpo_f_ice");         v2D = real(ylmo%tpo%now%f_ice,         c_double)
      case("tpo_f_grnd_bmb");    v2D = real(ylmo%tpo%now%f_grnd_bmb,    c_double)
      case("tpo_f_grnd_pin");    v2D = real(ylmo%tpo%now%f_grnd_pin,    c_double)
      case("tpo_dist_margin");   v2D = real(ylmo%tpo%now%dist_margin,   c_double)
      case("tpo_dist_grline");   v2D = real(ylmo%tpo%now%dist_grline,   c_double)
      case("tpo_dHidt_dyn_n");   v2D = real(ylmo%tpo%now%dHidt_dyn_n,   c_double)
      case("tpo_H_ice_n");       v2D = real(ylmo%tpo%now%H_ice_n,       c_double)
      case("tpo_z_srf_n");       v2D = real(ylmo%tpo%now%z_srf_n,       c_double)
      case("tpo_lsf_n");         v2D = real(ylmo%tpo%now%lsf_n,         c_double)
      case("tpo_H_ice_dyn");     v2D = real(ylmo%tpo%now%H_ice_dyn,     c_double)
      case("tpo_f_ice_dyn");     v2D = real(ylmo%tpo%now%f_ice_dyn,     c_double)
      case("tpo_tau_relax");     v2D = real(ylmo%tpo%now%tau_relax,     c_double)
      ! integer masks cast to double
      case("tpo_mask_adv");      v2D = real(ylmo%tpo%now%mask_adv,      c_double)
      case("tpo_mask_bed");      v2D = real(ylmo%tpo%now%mask_bed,      c_double)
      case("tpo_mask_grz");      v2D = real(ylmo%tpo%now%mask_grz,      c_double)
      case("tpo_mask_frnt");     v2D = real(ylmo%tpo%now%mask_frnt,     c_double)

      ! -----------------------------------------------------------------------
      ! ydyn%now — depth-integrated / surface / basal 2D fields
      ! -----------------------------------------------------------------------
      case("dyn_ux_bar");        v2D = real(ylmo%dyn%now%ux_bar,        c_double)
      case("dyn_uy_bar");        v2D = real(ylmo%dyn%now%uy_bar,        c_double)
      case("dyn_uxy_bar");       v2D = real(ylmo%dyn%now%uxy_bar,       c_double)
      case("dyn_ux_bar_prev");   v2D = real(ylmo%dyn%now%ux_bar_prev,   c_double)
      case("dyn_uy_bar_prev");   v2D = real(ylmo%dyn%now%uy_bar_prev,   c_double)
      case("dyn_ux_b");          v2D = real(ylmo%dyn%now%ux_b,          c_double)
      case("dyn_uy_b");          v2D = real(ylmo%dyn%now%uy_b,          c_double)
      case("dyn_uz_b");          v2D = real(ylmo%dyn%now%uz_b,          c_double)
      case("dyn_uxy_b");         v2D = real(ylmo%dyn%now%uxy_b,         c_double)
      case("dyn_ux_s");          v2D = real(ylmo%dyn%now%ux_s,          c_double)
      case("dyn_uy_s");          v2D = real(ylmo%dyn%now%uy_s,          c_double)
      case("dyn_uz_s");          v2D = real(ylmo%dyn%now%uz_s,          c_double)
      case("dyn_uxy_s");         v2D = real(ylmo%dyn%now%uxy_s,         c_double)
      case("dyn_ux_i_bar");      v2D = real(ylmo%dyn%now%ux_i_bar,      c_double)
      case("dyn_uy_i_bar");      v2D = real(ylmo%dyn%now%uy_i_bar,      c_double)
      case("dyn_uxy_i_bar");     v2D = real(ylmo%dyn%now%uxy_i_bar,     c_double)
      case("dyn_duxydt");        v2D = real(ylmo%dyn%now%duxydt,        c_double)
      case("dyn_duxdz_bar");     v2D = real(ylmo%dyn%now%duxdz_bar,     c_double)
      case("dyn_duydz_bar");     v2D = real(ylmo%dyn%now%duydz_bar,     c_double)
      case("dyn_taud_acx");      v2D = real(ylmo%dyn%now%taud_acx,      c_double)
      case("dyn_taud_acy");      v2D = real(ylmo%dyn%now%taud_acy,      c_double)
      case("dyn_taud");          v2D = real(ylmo%dyn%now%taud,          c_double)
      case("dyn_taub_acx");      v2D = real(ylmo%dyn%now%taub_acx,      c_double)
      case("dyn_taub_acy");      v2D = real(ylmo%dyn%now%taub_acy,      c_double)
      case("dyn_taub");          v2D = real(ylmo%dyn%now%taub,          c_double)
      case("dyn_taul_int_acx");  v2D = real(ylmo%dyn%now%taul_int_acx,  c_double)
      case("dyn_taul_int_acy");  v2D = real(ylmo%dyn%now%taul_int_acy,  c_double)
      case("dyn_qq_gl_acx");     v2D = real(ylmo%dyn%now%qq_gl_acx,     c_double)
      case("dyn_qq_gl_acy");     v2D = real(ylmo%dyn%now%qq_gl_acy,     c_double)
      case("dyn_qq_acx");        v2D = real(ylmo%dyn%now%qq_acx,        c_double)
      case("dyn_qq_acy");        v2D = real(ylmo%dyn%now%qq_acy,        c_double)
      case("dyn_qq");            v2D = real(ylmo%dyn%now%qq,            c_double)
      case("dyn_visc_eff_int");  v2D = real(ylmo%dyn%now%visc_eff_int,  c_double)
      case("dyn_N_eff");         v2D = real(ylmo%dyn%now%N_eff,         c_double)
      case("dyn_cb_tgt");        v2D = real(ylmo%dyn%now%cb_tgt,        c_double)
      case("dyn_cb_ref");        v2D = real(ylmo%dyn%now%cb_ref,        c_double)
      case("dyn_c_bed");         v2D = real(ylmo%dyn%now%c_bed,         c_double)
      case("dyn_beta_acx");      v2D = real(ylmo%dyn%now%beta_acx,      c_double)
      case("dyn_beta_acy");      v2D = real(ylmo%dyn%now%beta_acy,      c_double)
      case("dyn_beta");          v2D = real(ylmo%dyn%now%beta,          c_double)
      case("dyn_beta_eff");      v2D = real(ylmo%dyn%now%beta_eff,      c_double)
      case("dyn_f_vbvs");        v2D = real(ylmo%dyn%now%f_vbvs,        c_double)
      case("dyn_ssa_err_acx");   v2D = real(ylmo%dyn%now%ssa_err_acx,   c_double)
      case("dyn_ssa_err_acy");   v2D = real(ylmo%dyn%now%ssa_err_acy,   c_double)
      ! integer SSA masks cast to double
      case("dyn_ssa_mask_acx");  v2D = real(ylmo%dyn%now%ssa_mask_acx,  c_double)
      case("dyn_ssa_mask_acy");  v2D = real(ylmo%dyn%now%ssa_mask_acy,  c_double)
      ! strain_2D_class (dyn)
      case("dyn_strn2D_dxx");    v2D = real(ylmo%dyn%now%strn2D%dxx,    c_double)
      case("dyn_strn2D_dyy");    v2D = real(ylmo%dyn%now%strn2D%dyy,    c_double)
      case("dyn_strn2D_dxy");    v2D = real(ylmo%dyn%now%strn2D%dxy,    c_double)
      case("dyn_strn2D_dxz");    v2D = real(ylmo%dyn%now%strn2D%dxz,    c_double)
      case("dyn_strn2D_dyz");    v2D = real(ylmo%dyn%now%strn2D%dyz,    c_double)
      case("dyn_strn2D_de");     v2D = real(ylmo%dyn%now%strn2D%de,     c_double)
      case("dyn_strn2D_div");    v2D = real(ylmo%dyn%now%strn2D%div,    c_double)
      case("dyn_strn2D_f_shear"); v2D = real(ylmo%dyn%now%strn2D%f_shear, c_double)
      case("dyn_strn2D_eps_eig_1"); v2D = real(ylmo%dyn%now%strn2D%eps_eig_1, c_double)
      case("dyn_strn2D_eps_eig_2"); v2D = real(ylmo%dyn%now%strn2D%eps_eig_2, c_double)

      ! -----------------------------------------------------------------------
      ! ymat%now — depth-integrated 2D fields
      ! -----------------------------------------------------------------------
      case("mat_enh_bar");       v2D = real(ylmo%mat%now%enh_bar,       c_double)
      case("mat_ATT_bar");       v2D = real(ylmo%mat%now%ATT_bar,       c_double)
      case("mat_visc_bar");      v2D = real(ylmo%mat%now%visc_bar,      c_double)
      case("mat_visc_int");      v2D = real(ylmo%mat%now%visc_int,      c_double)
      case("mat_f_shear_bar");   v2D = real(ylmo%mat%now%f_shear_bar,   c_double)
      ! strain_2D_class (mat)
      case("mat_strn2D_dxx");    v2D = real(ylmo%mat%now%strn2D%dxx,    c_double)
      case("mat_strn2D_dyy");    v2D = real(ylmo%mat%now%strn2D%dyy,    c_double)
      case("mat_strn2D_dxy");    v2D = real(ylmo%mat%now%strn2D%dxy,    c_double)
      case("mat_strn2D_dxz");    v2D = real(ylmo%mat%now%strn2D%dxz,    c_double)
      case("mat_strn2D_dyz");    v2D = real(ylmo%mat%now%strn2D%dyz,    c_double)
      case("mat_strn2D_de");     v2D = real(ylmo%mat%now%strn2D%de,     c_double)
      case("mat_strn2D_div");    v2D = real(ylmo%mat%now%strn2D%div,    c_double)
      case("mat_strn2D_f_shear"); v2D = real(ylmo%mat%now%strn2D%f_shear, c_double)
      case("mat_strn2D_eps_eig_1"); v2D = real(ylmo%mat%now%strn2D%eps_eig_1, c_double)
      case("mat_strn2D_eps_eig_2"); v2D = real(ylmo%mat%now%strn2D%eps_eig_2, c_double)
      ! stress_2D_class (mat)
      case("mat_strs2D_txx");    v2D = real(ylmo%mat%now%strs2D%txx,    c_double)
      case("mat_strs2D_tyy");    v2D = real(ylmo%mat%now%strs2D%tyy,    c_double)
      case("mat_strs2D_txy");    v2D = real(ylmo%mat%now%strs2D%txy,    c_double)
      case("mat_strs2D_txz");    v2D = real(ylmo%mat%now%strs2D%txz,    c_double)
      case("mat_strs2D_tyz");    v2D = real(ylmo%mat%now%strs2D%tyz,    c_double)
      case("mat_strs2D_te");     v2D = real(ylmo%mat%now%strs2D%te,     c_double)
      case("mat_strs2D_tau_eig_1"); v2D = real(ylmo%mat%now%strs2D%tau_eig_1, c_double)
      case("mat_strs2D_tau_eig_2"); v2D = real(ylmo%mat%now%strs2D%tau_eig_2, c_double)

      ! -----------------------------------------------------------------------
      ! ytherm%now — 2D fields
      ! -----------------------------------------------------------------------
      case("thrm_f_pmp");        v2D = real(ylmo%thrm%now%f_pmp,        c_double)
      case("thrm_bmb_grnd");     v2D = real(ylmo%thrm%now%bmb_grnd,     c_double)
      case("thrm_Q_b");          v2D = real(ylmo%thrm%now%Q_b,          c_double)
      case("thrm_Q_ice_b");      v2D = real(ylmo%thrm%now%Q_ice_b,      c_double)
      case("thrm_T_prime_b");    v2D = real(ylmo%thrm%now%T_prime_b,    c_double)
      case("thrm_H_w");          v2D = real(ylmo%thrm%now%H_w,          c_double)
      case("thrm_dHwdt");        v2D = real(ylmo%thrm%now%dHwdt,        c_double)
      case("thrm_H_cts");        v2D = real(ylmo%thrm%now%H_cts,        c_double)
      case("thrm_Q_rock");       v2D = real(ylmo%thrm%now%Q_rock,       c_double)

      case DEFAULT
        write(*,*) "yelmo_get_var2D:: variable not found: "//trim(f_name)// &
                 " — returning MISSING_VALUE"
        v2D = real(MISSING_VALUE_DEFAULT, c_double)
    end select

    nullify(ylmo)

  end subroutine yelmo_get_var2D


  ! =========================================================================
  ! 3D variable getter  (nx, ny, nz — caller must know the correct nz)
  ! =========================================================================
  subroutine yelmo_get_var3D(v3D, nx, ny, nz, name, alias) bind(C, name="yelmo_get_var3D")
    use iso_c_binding
    integer(c_int), value         :: nx, ny, nz
    character(c_char), intent(in) :: name(*)
    real(c_double), intent(out)   :: v3D(nx, ny, nz)
    character(c_char), intent(in) :: alias(*)

    ! Local variables
    character(len=56) :: f_name

    ! Convert C string to Fortran string (null-terminated scan)
    f_name = trim(c_to_f_string(name))

    call yelmo_set_alias(alias)

    select case(trim(f_name))

      ! -----------------------------------------------------------------------
      ! ydyn%now — full 3D velocity / strain fields
      ! -----------------------------------------------------------------------
      case("dyn_ux");            v3D = real(ylmo%dyn%now%ux,            c_double)
      case("dyn_uy");            v3D = real(ylmo%dyn%now%uy,            c_double)
      case("dyn_uxy");           v3D = real(ylmo%dyn%now%uxy,           c_double)
      case("dyn_uz");            v3D = real(ylmo%dyn%now%uz,            c_double)
      case("dyn_uz_star");       v3D = real(ylmo%dyn%now%uz_star,       c_double)
      case("dyn_ux_i");          v3D = real(ylmo%dyn%now%ux_i,          c_double)
      case("dyn_uy_i");          v3D = real(ylmo%dyn%now%uy_i,          c_double)
      case("dyn_duxdz");         v3D = real(ylmo%dyn%now%duxdz,         c_double)
      case("dyn_duydz");         v3D = real(ylmo%dyn%now%duydz,         c_double)
      case("dyn_de_eff");        v3D = real(ylmo%dyn%now%de_eff,        c_double)
      case("dyn_visc_eff");      v3D = real(ylmo%dyn%now%visc_eff,      c_double)
      ! strain_3D_class (dyn)
      case("dyn_strn_dxx");      v3D = real(ylmo%dyn%now%strn%dxx,      c_double)
      case("dyn_strn_dyy");      v3D = real(ylmo%dyn%now%strn%dyy,      c_double)
      case("dyn_strn_dxy");      v3D = real(ylmo%dyn%now%strn%dxy,      c_double)
      case("dyn_strn_dxz");      v3D = real(ylmo%dyn%now%strn%dxz,      c_double)
      case("dyn_strn_dyz");      v3D = real(ylmo%dyn%now%strn%dyz,      c_double)
      case("dyn_strn_de");       v3D = real(ylmo%dyn%now%strn%de,       c_double)
      case("dyn_strn_div");      v3D = real(ylmo%dyn%now%strn%div,      c_double)
      case("dyn_strn_f_shear");  v3D = real(ylmo%dyn%now%strn%f_shear,  c_double)
      ! jacobian_3D_class (dyn)
      case("dyn_jvel_dxx");      v3D = real(ylmo%dyn%now%jvel%dxx,      c_double)
      case("dyn_jvel_dxy");      v3D = real(ylmo%dyn%now%jvel%dxy,      c_double)
      case("dyn_jvel_dxz");      v3D = real(ylmo%dyn%now%jvel%dxz,      c_double)
      case("dyn_jvel_dyx");      v3D = real(ylmo%dyn%now%jvel%dyx,      c_double)
      case("dyn_jvel_dyy");      v3D = real(ylmo%dyn%now%jvel%dyy,      c_double)
      case("dyn_jvel_dyz");      v3D = real(ylmo%dyn%now%jvel%dyz,      c_double)
      case("dyn_jvel_dzx");      v3D = real(ylmo%dyn%now%jvel%dzx,      c_double)
      case("dyn_jvel_dzy");      v3D = real(ylmo%dyn%now%jvel%dzy,      c_double)
      case("dyn_jvel_dzz");      v3D = real(ylmo%dyn%now%jvel%dzz,      c_double)

      ! -----------------------------------------------------------------------
      ! ymat%now — 3D material fields
      ! -----------------------------------------------------------------------
      case("mat_enh");           v3D = real(ylmo%mat%now%enh,           c_double)
      case("mat_enh_bnd");       v3D = real(ylmo%mat%now%enh_bnd,       c_double)
      case("mat_ATT");           v3D = real(ylmo%mat%now%ATT,           c_double)
      case("mat_visc");          v3D = real(ylmo%mat%now%visc,          c_double)
      case("mat_dep_time");      v3D = real(ylmo%mat%now%dep_time,      c_double)
      case("mat_depth_iso");     v3D = real(ylmo%mat%now%depth_iso,     c_double)
      ! strain_3D_class (mat)
      case("mat_strn_dxx");      v3D = real(ylmo%mat%now%strn%dxx,      c_double)
      case("mat_strn_dyy");      v3D = real(ylmo%mat%now%strn%dyy,      c_double)
      case("mat_strn_dxy");      v3D = real(ylmo%mat%now%strn%dxy,      c_double)
      case("mat_strn_dxz");      v3D = real(ylmo%mat%now%strn%dxz,      c_double)
      case("mat_strn_dyz");      v3D = real(ylmo%mat%now%strn%dyz,      c_double)
      case("mat_strn_de");       v3D = real(ylmo%mat%now%strn%de,       c_double)
      case("mat_strn_div");      v3D = real(ylmo%mat%now%strn%div,      c_double)
      case("mat_strn_f_shear");  v3D = real(ylmo%mat%now%strn%f_shear,  c_double)
      ! stress_3D_class (mat)
      case("mat_strs_txx");      v3D = real(ylmo%mat%now%strs%txx,      c_double)
      case("mat_strs_tyy");      v3D = real(ylmo%mat%now%strs%tyy,      c_double)
      case("mat_strs_txy");      v3D = real(ylmo%mat%now%strs%txy,      c_double)
      case("mat_strs_txz");      v3D = real(ylmo%mat%now%strs%txz,      c_double)
      case("mat_strs_tyz");      v3D = real(ylmo%mat%now%strs%tyz,      c_double)
      case("mat_strs_te");       v3D = real(ylmo%mat%now%strs%te,       c_double)

      ! -----------------------------------------------------------------------
      ! ytherm%now — 3D thermodynamic fields
      ! -----------------------------------------------------------------------
      case("thrm_enth");         v3D = real(ylmo%thrm%now%enth,         c_double)
      case("thrm_T_ice");        v3D = real(ylmo%thrm%now%T_ice,        c_double)
      case("thrm_omega");        v3D = real(ylmo%thrm%now%omega,        c_double)
      case("thrm_T_pmp");        v3D = real(ylmo%thrm%now%T_pmp,        c_double)
      case("thrm_T_prime");      v3D = real(ylmo%thrm%now%T_prime,      c_double)
      case("thrm_Q_strn");       v3D = real(ylmo%thrm%now%Q_strn,       c_double)
      case("thrm_dQsdT");        v3D = real(ylmo%thrm%now%dQsdT,        c_double)
      case("thrm_cp");           v3D = real(ylmo%thrm%now%cp,           c_double)
      case("thrm_kt");           v3D = real(ylmo%thrm%now%kt,           c_double)
      case("thrm_advecxy");      v3D = real(ylmo%thrm%now%advecxy,      c_double)
      case("thrm_enth_rock");    v3D = real(ylmo%thrm%now%enth_rock,    c_double)
      case("thrm_T_rock");       v3D = real(ylmo%thrm%now%T_rock,       c_double)

      case DEFAULT
        write(*,*) "yelmo_get_var3D:: variable not found: "//trim(f_name)// &
                 " — returning MISSING_VALUE"
        v3D = real(MISSING_VALUE_DEFAULT, c_double)
    end select

    nullify(ylmo)

  end subroutine yelmo_get_var3D

! =========================================================================
  ! C-bound boundary field setter: restricted to ybound_class fields only.
  ! =========================================================================
  subroutine yelmo_set_var2D(v2D, nx, ny, name, alias) bind(C, name="yelmo_set_var2D")
    use iso_c_binding
    integer(c_int),    value      :: nx, ny
    real(c_double),    intent(in) :: v2D(nx, ny)
    character(c_char), intent(in) :: name(*)
    character(c_char), intent(in) :: alias(*)

    ! Local variables
    character(len=256) :: f_name

    f_name = trim(c_to_f_string(name))

    call yelmo_set_alias(alias)
    
    select case(trim(f_name))
      case("bnd_z_bed");       ylmo%bnd%z_bed      = real(v2D, wp)
      case("bnd_z_bed_sd");    ylmo%bnd%z_bed_sd   = real(v2D, wp)
      case("bnd_z_sl");        ylmo%bnd%z_sl       = real(v2D, wp)
      case("bnd_H_sed");       ylmo%bnd%H_sed      = real(v2D, wp)
      case("bnd_smb");         ylmo%bnd%smb        = real(v2D, wp)
      case("bnd_T_srf");       ylmo%bnd%T_srf      = real(v2D, wp)
      case("bnd_bmb_shlf");    ylmo%bnd%bmb_shlf   = real(v2D, wp)
      case("bnd_fmb_shlf");    ylmo%bnd%fmb_shlf   = real(v2D, wp)
      case("bnd_T_shlf");      ylmo%bnd%T_shlf     = real(v2D, wp)
      case("bnd_Q_geo");       ylmo%bnd%Q_geo      = real(v2D, wp)
      case("bnd_enh_srf");     ylmo%bnd%enh_srf    = real(v2D, wp)
      case("bnd_H_ice_ref");   ylmo%bnd%H_ice_ref  = real(v2D, wp)
      case("bnd_z_bed_ref");   ylmo%bnd%z_bed_ref  = real(v2D, wp)
      case("bnd_tau_relax");   ylmo%bnd%tau_relax  = real(v2D, wp)

      case DEFAULT
        write(*,*) "yelmo_set_var2D:: variable not found or not settable: "//trim(f_name)
    end select

    nullify(ylmo)

  end subroutine yelmo_set_var2D

  function c_to_f_string(c_str) result(f_str)
    use iso_c_binding
    character(c_char), intent(in) :: c_str(*)
    character(len=1028)           :: f_str
    integer                       :: i

    f_str = " "
    do i = 1, 1028
      if (c_str(i) == c_null_char) exit
      f_str(i:i) = c_str(i)
    end do
  end function
  
end module yelmo_c_api