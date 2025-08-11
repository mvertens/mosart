module mosart_histflds

   ! Module containing initialization of history fields and files
   ! This is the module that the user must modify in order to add new
   ! history fields or modify defaults associated with existing history
   ! fields.

   use shr_kind_mod    , only : r8 => shr_kind_r8
   use mosart_histfile , only : mosart_hist_addfld, mosart_hist_printflds
   use mosart_data     , only : ctl, Trunoff

   implicit none
   private

   public :: mosart_histflds_init
   public :: mosart_histflds_set

  type, public ::  hist_pointer_type
     real(r8), pointer :: data(:) => null()
  end type hist_pointer_type

  type(hist_pointer_type), allocatable :: h_runofflnd(:)
  type(hist_pointer_type), allocatable :: h_runoffocn(:)
  type(hist_pointer_type), allocatable :: h_runofftot(:)
  type(hist_pointer_type), allocatable :: h_direct(:)
  type(hist_pointer_type), allocatable :: h_dvolrdtlnd(:)
  type(hist_pointer_type), allocatable :: h_dvolrdtocn(:)
  type(hist_pointer_type), allocatable :: h_volr(:)
  type(hist_pointer_type), allocatable :: h_qsur_liq_nonh2o(:)

  type(hist_pointer_type) :: h_qsur_liq
  type(hist_pointer_type) :: h_qsur_ice
  type(hist_pointer_type) :: h_qsub_liq
  type(hist_pointer_type) :: h_qgwl_liq
  type(hist_pointer_type) :: h_direct_glc_liq
  type(hist_pointer_type) :: h_direct_glc_ice
  type(hist_pointer_type) :: h_volr_mch
  type(hist_pointer_type) :: h_qglc_liq_input
  type(hist_pointer_type) :: h_qglc_ice_input

!------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

   subroutine mosart_histflds_init()

      ! Arguments

      ! Local variables
      integer :: nt
      integer :: begr
      integer :: endr
      integer :: ntracers_tot
      integer :: ntracers_nonh2o

      begr = ctl%begr
      endr = ctl%endr
      ntracers_tot = ctl%ntracers_tot
      ntracers_nonh2o = ctl%ntracers_nonh2o

      !-------------------------------------------------------
      ! Allocate memory for module variables
      !-------------------------------------------------------

      ! Output

      allocate(h_runofflnd(ntracers_tot))
      allocate(h_runoffocn(ntracers_tot))
      allocate(h_runofftot(ntracers_tot))
      allocate(h_direct(ntracers_tot))
      allocate(h_dvolrdtlnd(ntracers_tot))
      allocate(h_dvolrdtocn(ntracers_tot))
      allocate(h_volr(ntracers_tot))

      do nt = 1,ntracers_tot
         allocate(h_runofflnd(nt)%data(begr:endr))
         allocate(h_runoffocn(nt)%data(begr:endr))
         allocate(h_runofftot(nt)%data(begr:endr))
         allocate(h_direct(nt)%data(begr:endr))
         allocate(h_volr(nt)%data(begr:endr))
         allocate(h_dvolrdtlnd(nt)%data(begr:endr))
         allocate(h_dvolrdtocn(nt)%data(begr:endr))
      end do
      allocate(h_volr_mch%data(begr:endr))
      allocate(h_direct_glc_liq%data(begr:endr))
      allocate(h_direct_glc_ice%data(begr:endr))

      ! Input

      allocate(h_qsur_liq%data(begr:endr))
      allocate(h_qsur_ice%data(begr:endr))
      allocate(h_qsub_liq%data(begr:endr))
      allocate(h_qgwl_liq%data(begr:endr))
      allocate(h_qglc_liq_input%data(begr:endr))
      allocate(h_qglc_ice_input%data(begr:endr))

      allocate(h_qsur_liq_nonh2o(ntracers_nonh2o))
      do nt = 1,ctl%ntracers_nonh2o
         allocate(h_qsur_liq_nonh2o(nt)%data(begr:endr))
      end do

      !-------------------------------------------------------
      ! Build master field list of all possible fields in a history file.
      ! Each field has associated with it a ``long\_name'' netcdf attribute that
      ! describes what the field is, and a ``units'' attribute. A subroutine is
      ! called to add each field to the masterlist.
      !-------------------------------------------------------

      ! ---------------
      ! Output
      ! ---------------

      do nt = 1,ctl%ntracers_tot
         call mosart_hist_addfld (fname='RIVER_DISCHARGE_OVER_LAND'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART river basin flow: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_runofflnd(nt)%data, default='active')

         call mosart_hist_addfld (fname='RIVER_DISCHARGE_TO_OCEAN'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART river discharge into ocean: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_runoffocn(nt)%data, default='active')

         call mosart_hist_addfld (fname='DIRECT_DISCHARGE_TO_OCEAN'//'_'//trim(ctl%tracer_names(nt)), units='m3/s', &
              avgflag='A', long_name='MOSART direct discharge into ocean: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_direct(nt)%data, default='active')

         call mosart_hist_addfld (fname='TOTAL_DISCHARGE_TO_OCEAN'//'_'//trim(ctl%tracer_names(nt)), units='m3/s', &
              avgflag='A', long_name='MOSART total discharge into ocean: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_runofftot(nt)%data, default='active')

         call mosart_hist_addfld (fname='STORAGE'//'_'//trim(ctl%tracer_names(nt)), units='m3',  &
              avgflag='A', long_name='MOSART storage: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_volr(nt)%data, default='active')

         call mosart_hist_addfld (fname='DVOLRDT_LND'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART land change in storage: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_dvolrdtlnd(nt)%data, default='active')

         call mosart_hist_addfld (fname='DVOLRDT_OCN'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART ocean change of storage: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_dvolrdtocn(nt)%data, default='active')

      end do

      call mosart_hist_addfld (fname='DIRECT_DISCHARGE_TO_OCEAN_GLC_LIQ', units='m3/s', &
           avgflag='A', long_name='MOSART direct discharge into ocean from glc liq: ', &
           ptr_rof=h_direct_glc_liq%data, default='active')

      call mosart_hist_addfld (fname='DIRECT_DISCHARGE_TO_OCEAN_GLC_ICE', units='m3/s', &
           avgflag='A', long_name='MOSART direct discharge into ocean from glc ice: ', &
           ptr_rof=h_direct_glc_ice%data, default='active')

      ! MOSART (unlike the CLM) does not have the history_tape_in_use
      ! capability, so both models throw an error when h0i is empty. For this
      ! reason MOSART always need at least one instantaneous field so
      ! that h0i will not be empty.
      call mosart_hist_addfld (fname='STORAGE_MCH', units='m3',  &
           avgflag='I', long_name='MOSART main channelstorage', &
           ptr_rof=h_volr_mch%data, default='active')

      call mosart_hist_addfld (fname='QIRRIG_FROM_COUPLER', units='m3/s',  &
           avgflag='A', long_name='Amount of water used for irrigation (total flux received from coupler)', &
           ptr_rof=ctl%qirrig_liq, default='active')

      call mosart_hist_addfld (fname='QIRRIG_ACTUAL', units='m3/s',  &
           avgflag='A', long_name='Actual irrigation (if limited by river storage)', &
           ptr_rof=ctl%qirrig_actual, default='active')

      ! ---------------
      ! Input
      ! ---------------

      call mosart_hist_addfld (fname='QSUR_LIQ_INPUT', units='m3/s',  &
           avgflag='A', long_name='MOSART input surface runoff liquid water from land', &
           ptr_rof=h_qsur_liq%data, default='active')

      call mosart_hist_addfld (fname='QSUR_ICE_INPUT', units='m3/s',  &
           avgflag='A', long_name='MOSART input surface runoff ice from land', &
           ptr_rof=h_qsur_ice%data, default='active')

      do nt = 1,ctl%ntracers_nonh2o
         call mosart_hist_addfld (fname='QSUR_NONH2O'//'_'//trim(ctl%tracer_names(nt+2))//'_INPUT', units='m3/s',  &
              avgflag='A', long_name='MOSART input surface runoff non-h2o tracers from land', &
              ptr_rof=h_qsur_liq_nonh2o(nt)%data, default='active')
      end do

      call mosart_hist_addfld (fname='QSUB_LIQ_INPUT', units='m3/s',  &
           avgflag='A', long_name='MOSART input subsurface liquid water runoff from land', &
           ptr_rof=h_qsub_liq%data, default='active')

      call mosart_hist_addfld (fname='QGWL_LIQ_INPUT', units='m3/s',  &
           avgflag='A', long_name='MOSART input glacier/wetland/lake runoff from land', &
           ptr_rof=h_qgwl_liq%data, default='active')

      call mosart_hist_addfld (fname='QGLC_LIQ_INPUT', units='m3',  &
           avgflag='A', long_name='MOSART input liquid runoff from glacier', &
           ptr_rof=h_qglc_liq_input%data, default='active')

      call mosart_hist_addfld (fname='QGLC_ICE_INPUT', units='m3',  &
           avgflag='A', long_name='MOSART input ice runoff from glacier', &
           ptr_rof=h_qglc_ice_input%data, default='active')

      ! print masterlist of history fields
      call mosart_hist_printflds()

   end subroutine mosart_histflds_init

   !-----------------------------------------------------------------------

   subroutine mosart_histflds_set()

      !-----------------------------------------------------------------------
      ! Set mosart history fields as 1d pointer arrays
      !-----------------------------------------------------------------------

      ! Local variables
      integer :: nt
      integer :: nt_liq, nt_ice
      integer :: ntracers_tot, ntracers_nonh2o

      ntracers_tot = ctl%ntracers_tot
      ntracers_nonh2o = ctl%ntracers_nonh2o
      nt_liq = ctl%nt_liq
      nt_ice = ctl%nt_ice

      ! Input
      h_qsur_liq%data(:) = ctl%qsur_liq(:)
      h_qsur_ice%data(:) = ctl%qsur_ice(:)
      h_qsub_liq%data(:) = ctl%qsub_liq(:)
      h_qgwl_liq%data(:) = ctl%qgwl_liq(:)
      h_qglc_liq_input%data(:) = ctl%qglc_liq(:)
      h_qglc_ice_input%data(:) = ctl%qglc_ice(:)
      do nt = 1,ctl%ntracers_nonh2o
         h_qsur_liq_nonh2o(nt)%data(:) = ctl%qsur_liq_nonh2o(:,nt)
      end do

      ! Output
      do nt = 1,ntracers_tot
         h_runofflnd(nt)%data(:)  = ctl%runofflnd(:,nt)
         h_runoffocn(nt)%data(:)  = ctl%runoffocn(:,nt)
         h_runofftot(nt)%data(:)  = ctl%runofftot(:,nt)
         h_direct(nt)%data(:)     = ctl%direct(:,nt)
         h_volr(nt)%data(:)       = ctl%volr(:,nt)
         h_dvolrdtlnd(nt)%data(:) = ctl%dvolrdtlnd(:,nt)
         h_dvolrdtocn(nt)%data(:) = ctl%dvolrdtocn(:,nt)
      end do
      h_volr_mch%data(:) = Trunoff%wr(:,1)
      h_direct_glc_liq%data(:) = ctl%direct_glc(:,nt_liq)
      h_direct_glc_ice%data(:) = ctl%direct_glc(:,nt_ice)

   end subroutine mosart_histflds_set

end module mosart_histflds
