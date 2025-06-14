module rof_import_export

  use ESMF                , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF                , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF                , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF                , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF                , only : ESMF_StateGet, ESMF_FieldGet, ESMF_Field
  use ESMF                , only : operator(/=), operator(==)
  use NUOPC               , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model         , only : NUOPC_ModelGet
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_abort
  use shr_infnan_mod      , only : isnan => shr_infnan_isnan
  use mosart_vars         , only : iulog, mainproc, mpicom_rof, ice_runoff
  use mosart_data         , only : ctl, TRunoff, TUnit
  use mosart_timemanager  , only : get_nstep
  use nuopc_shr_methods   , only : chkerr

  implicit none
  private ! except

  public  :: advertise_fields
  public  :: realize_fields
  public  :: import_fields
  public  :: export_fields

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_getimport1d
  private :: state_getimport2d
  private :: state_setexport1d
  private :: state_setexport2d
  private :: check_for_nans
  private :: fldchk

  interface check_for_nans
     module procedure check_for_nans_1d
     module procedure check_for_nans_2d
  end interface check_for_nans

  type fld_list_type
     character(len=128) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToRof_num = 0
  integer                :: fldsFrRof_num = 0
  logical                :: flds_r2l_stream_channel_depths = .false.   ! If should pass the channel depth fields needed for the hillslope model
  type (fld_list_type)   :: fldsToRof(fldsMax)
  type (fld_list_type)   :: fldsFrRof(fldsMax)

  ! area correction factors for fluxes send and received from mediator
  real(r8), allocatable :: mod2med_areacor(:)
  real(r8), allocatable :: med2mod_areacor(:)

  character(*),parameter :: F01 = "('(mosart_import_export) ',a,i5,2x,i8,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, flds_scalar_name, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)       :: importState
    type(ESMF_State)       :: exportState
    character(ESMF_MAXSTR) :: cvalue          ! Character string read from driver attribute
    logical                :: isPresent       ! Atribute is present
    logical                :: isSet           ! Atribute is set
    integer                :: n
    integer                :: ntracers_nonh2o
    character(len=*), parameter :: subname='(rof_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ntracers_nonh2o = ctl%ntracers_nonh2o

    !--------------------------------
    ! Advertise export fields
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="flds_r2l_stream_channel_depths", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) read(cvalue,*) flds_r2l_stream_channel_depths

    call fldlist_add(fldsFrRof_num, fldsFrRof, trim(flds_scalar_name))
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofl')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofi')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofl_glc')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofi_glc')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_flood')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_volr')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_volrmch')
    if ( flds_r2l_stream_channel_depths )then
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Sr_tdepth')
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Sr_tdepth_max')
    end if
    if (ntracers_nonh2o > 1) then
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofl_nonh2o', ungridded_lbound=1, ungridded_ubound=ntracers_nonh2o)
    else if (ntracers_nonh2o == 1) then
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofl_nonh2o')
    end if

    do n = 1,fldsFrRof_num
       call NUOPC_Advertise(exportState, standardName=fldsFrRof(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    call fldlist_add(fldsToRof_num, fldsToRof, trim(flds_scalar_name))
    if (ntracers_nonh2o > 1) then
       call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsur_nonh2o', ungridded_lbound=1, ungridded_ubound=ntracers_nonh2o)
    else  if (ntracers_nonh2o == 1) then
       call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsur_nonh2o')
    end if
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsur')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofgwl')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsub')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofi')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_irrig')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Fgrg_rofl') ! liq runoff from glc
    call fldlist_add(fldsToRof_num, fldsToRof, 'Fgrg_rofi') ! ice runoff from glc

    do n = 1,fldsToRof_num
       call NUOPC_Advertise(importState, standardName=fldsToRof(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine advertise_fields

  !===============================================================================
  subroutine realize_fields(gcomp, Emesh, flds_scalar_name, flds_scalar_num, rc)

    use ESMF          , only : ESMF_GridComp, ESMF_StateGet
    use ESMF          , only : ESMF_Mesh, ESMF_MeshGet
    use ESMF          , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldRegridGetArea
    use shr_const_mod , only : shr_const_rearth
    use shr_mpi_mod   , only : shr_mpi_min, shr_mpi_max

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_Mesh)     , intent(in)    :: Emesh
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_State)      :: importState
    type(ESMF_State)      :: exportState
    type(ESMF_Field)      :: lfield
    integer               :: numOwnedElements
    integer               :: n,g
    real(r8), allocatable :: mesh_areas(:)
    real(r8), allocatable :: model_areas(:)
    real(r8), pointer     :: dataptr(:)
    real(r8)              :: re = shr_const_rearth*0.001_r8 ! radius of earth (km)
    real(r8)              :: max_mod2med_areacor
    real(r8)              :: max_med2mod_areacor
    real(r8)              :: min_mod2med_areacor
    real(r8)              :: min_med2mod_areacor
    real(r8)              :: max_mod2med_areacor_glob
    real(r8)              :: max_med2mod_areacor_glob
    real(r8)              :: min_mod2med_areacor_glob
    real(r8)              :: min_med2mod_areacor_glob
    character(len=*), parameter :: subname='(rof_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrRof, &
         numflds=fldsFrRof_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':MosartExport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToRof, &
         numflds=fldsToRof_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':MosartImport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine areas for regridding
    call ESMF_MeshGet(Emesh, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_StateGet(exportState, itemName=trim(fldsFrRof(2)%stdname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(mesh_areas(numOwnedElements))
    mesh_areas(:) = dataptr(:)

    ! Determine model areas
    allocate(model_areas(numOwnedElements))
    allocate(mod2med_areacor(numOwnedElements))
    allocate(med2mod_areacor(numOwnedElements))
    n = 0
    do g = ctl%begr,ctl%endr
       n = n + 1
       model_areas(n) = ctl%area(g)*1.0e-6_r8/(re*re)
       mod2med_areacor(n) = model_areas(n) / mesh_areas(n)
       med2mod_areacor(n) = mesh_areas(n) / model_areas(n)
    end do
    deallocate(model_areas)
    deallocate(mesh_areas)

    min_mod2med_areacor = minval(mod2med_areacor)
    max_mod2med_areacor = maxval(mod2med_areacor)
    min_med2mod_areacor = minval(med2mod_areacor)
    max_med2mod_areacor = maxval(med2mod_areacor)
    call shr_mpi_max(max_mod2med_areacor, max_mod2med_areacor_glob, mpicom_rof)
    call shr_mpi_min(min_mod2med_areacor, min_mod2med_areacor_glob, mpicom_rof)
    call shr_mpi_max(max_med2mod_areacor, max_med2mod_areacor_glob, mpicom_rof)
    call shr_mpi_min(min_med2mod_areacor, min_med2mod_areacor_glob, mpicom_rof)

    if (mainproc) then
       write(iulog,'(2A,2g23.15,A )') trim(subname),' :  min_mod2med_areacor, max_mod2med_areacor ',&
            min_mod2med_areacor_glob, max_mod2med_areacor_glob, 'MOSART'
       write(iulog,'(2A,2g23.15,A )') trim(subname),' :  min_med2mod_areacor, max_med2mod_areacor ',&
            min_med2mod_areacor_glob, max_med2mod_areacor_glob, 'MOSART'
    end if

    if (fldchk(importState, 'Fgrg_rofl') .and. fldchk(importState, 'Fgrg_rofl')) then
      ctl%rof_from_glc = .true.
    else
      ctl%rof_from_glc = .false.
    end if
    if (mainproc) then
      write(iulog,'(A,l1)') trim(subname) //' rof from glc is ',ctl%rof_from_glc
    end if

  end subroutine realize_fields

  !===============================================================================
  subroutine import_fields( gcomp, begr, endr, rc )

    !---------------------------------------------------------------------------
    ! Obtain the runoff input from the mediator and convert from kg/m2s to m3/s
    !---------------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(in)  :: begr, endr
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State) :: importState
    integer          :: ntracers_tot
    integer          :: ntracers_nonh2o
    character(len=*), parameter :: subname='(rof_import_export:import_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get import state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ntracers_tot = ctl%ntracers_tot
    ntracers_nonh2o = ctl%ntracers_nonh2o

    ! determine output array and scale by unit convertsion
    ! NOTE: the call to state_getimport will convert from input kg/m2/s to m3/s

    if (ntracers_nonh2o > 1) then
       call state_getimport2d(importState, 'Flrl_rofsur_nonh2o', begr, endr, ntracers_nonh2o, ctl%area, &
            output=ctl%qsur_liq_nonh2o, do_area_correction=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (ntracers_nonh2o == 1) then
       call state_getimport1d(importState, 'Flrl_rofsur_nonh2o', begr, endr, ctl%area, &
            output=ctl%qsur_liq_nonh2o(:,1), do_area_correction=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call state_getimport1d(importState, 'Flrl_rofsur', begr,endr, ctl%area, output=ctl%qsur_liq, &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport1d(importState, 'Flrl_rofsub', begr, endr, ctl%area, output=ctl%qsub_liq, &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport1d(importState, 'Flrl_rofgwl', begr, endr, ctl%area, output=ctl%qgwl_liq, &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport1d(importState, 'Flrl_rofi', begr, endr, ctl%area, output=ctl%qsur_ice(:), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport1d(importState, 'Flrl_irrig', begr, endr, ctl%area, output=ctl%qirrig_liq, &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ctl%rof_from_glc) then
      call state_getimport1d(importState, 'Fgrg_rofl', begr, endr, ctl%area, output=ctl%qglc_liq, &
           do_area_correction=.true., rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call state_getimport1d(importState, 'Fgrg_rofi', begr, endr, ctl%area, output=ctl%qglc_ice, &
           do_area_correction=.true., rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      ctl%qglc_liq(:) = 0._r8
      ctl%qglc_ice(:) = 0._r8
    end if

  end subroutine import_fields

  !====================================================================================
  subroutine export_fields (gcomp, begr, endr, ntracers_nonh2o, rc)

    !---------------------------------------------------------------------------
    ! Send the runoff model export state to the mediator and convert from m3/s to kg/m2s
    !---------------------------------------------------------------------------

    ! input/output/variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(in)  :: begr, endr
    integer, intent(in)  :: ntracers_nonh2o
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State) :: exportState
    integer          :: nr,nt
    integer          :: nliq, nice
    real(r8)         :: rof_liq(begr:endr)
    real(r8)         :: rof_ice(begr:endr)
    real(r8)         :: rof_liq_nonh2o(begr:endr,ntracers_nonh2o)
    real(r8)         :: rof_liq_glc(begr:endr)
    real(r8)         :: rof_ice_glc(begr:endr)
    real(r8)         :: flood(begr:endr)
    real(r8)         :: volr(begr:endr)
    real(r8)         :: volrmch(begr:endr)
    real(r8)         :: tdepth(begr:endr)
    real(r8)         :: tdepth_max(begr:endr)
    logical, save    :: first_time = .true.
    character(len=*), parameter :: subname='(rof_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (first_time) then
       if (mainproc) then
          if ( ice_runoff )then
             write(iulog,*)'Snow capping will flow out in frozen river runoff'
          else
             write(iulog,*)'Snow capping will flow out in liquid river runoff'
          endif
       endif
       first_time = .false.
    end if

    ! Set tracers
    nliq = ctl%nt_liq
    nice = ctl%nt_ice

    ! water tracers from lnd
    if ( ice_runoff )then
       ! separate water liquid and ice runoff
       do nr = begr,endr
          rof_liq(nr) = ctl%direct(nr,nliq) / (ctl%area(nr)*0.001_r8)
          rof_ice(nr) = ctl%direct(nr,nice) / (ctl%area(nr)*0.001_r8)
          if (ctl%mask(nr) >= 2) then
             ! water liquid and ice runoff are treated separately - this is what goes to the ocean
             rof_liq(nr) = rof_liq(nr) + ctl%runoff(nr,nliq) / (ctl%area(nr)*0.001_r8)
             rof_ice(nr) = rof_ice(nr) + ctl%runoff(nr,nice) / (ctl%area(nr)*0.001_r8)
          end if
       end do
    else
       ! ice runoff added to water liquid runoff, ice runoff is zero
       do nr = begr,endr
          rof_liq(nr) = (ctl%direct(nr,nice) + ctl%direct(nr,nliq)) / (ctl%area(nr)*0.001_r8)
          if (ctl%mask(nr) >= 2) then
             rof_liq(nr) = rof_liq(nr) + (ctl%runoff(nr,nice) + ctl%runoff(nr,nliq)) / (ctl%area(nr)*0.001_r8)
          endif
          rof_ice(nr) = 0._r8
       end do
    end if

    ! non-water liquid tracers from lnd to ocn
    do nt = 1,ntracers_nonh2o
       do nr = begr,endr
          if (ctl%mask(nr) >= 2) then
             rof_liq_nonh2o(nr,nt) = ctl%runoff(nr,nt+2) / (ctl%area(nr)*0.001_r8)
          else
             rof_liq_nonh2o(nr,nt) = 0._r8
          endif
       end do
    end do

    ! water tracers from glc
    do nr = begr,endr
      rof_liq_glc(nr) = ctl%direct_glc(nr,nliq) / (ctl%area(nr)*0.001_r8)
      rof_ice_glc(nr) = ctl%direct_glc(nr,nice) / (ctl%area(nr)*0.001_r8)
    end do

    ! Flooding back to land, sign convention is positive in land->rof direction
    ! so if water is sent from rof to land, the flux must be negative.
    ! scs: is there a reason for the wr+wt rather than volr (wr+wt+wh)?
    ! volr(n) = (Trunoff%wr(n,nliq) + Trunoff%wt(n,nliq)) / ctl%area(n)

    do nr = begr, endr
       flood(nr)   = -ctl%flood(nr) / (ctl%area(nr)*0.001_r8)
       volr(nr)    =  ctl%volr(nr,nliq) / ctl%area(nr)
       volrmch(nr) =  Trunoff%wr(nr,nliq) / ctl%area(nr)
       if ( flds_r2l_stream_channel_depths )then
          tdepth(nr)  = Trunoff%yt(nr,nliq)
          ! assume height to width ratio is the same for tributaries and main channel
          tdepth_max(nr) = max(TUnit%twidth0(nr),0._r8)*(TUnit%rdepth(nr)/TUnit%rwidth(nr))
        end if
    end do

    if (ntracers_nonh2o > 1) then
       call state_setexport2d(exportState, 'Forr_rofl_nonh2o', begr, endr, ntracers_nonh2o, input=rof_liq_nonh2o, &
            do_area_correction=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (ntracers_nonh2o == 1) then
       call state_setexport1d(exportState, 'Forr_rofl_nonh2o', begr, endr, input=rof_liq_nonh2o(:,1), &
            do_area_correction=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call state_setexport1d(exportState, 'Forr_rofl', begr, endr, input=rof_liq, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport1d(exportState, 'Forr_rofi', begr, endr, input=rof_ice, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport1d(exportState, 'Forr_rofl_glc', begr, endr, input=rof_liq_glc, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport1d(exportState, 'Forr_rofi_glc', begr, endr, input=rof_ice_glc, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport1d(exportState, 'Flrr_flood', begr, endr, input=flood, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport1d(exportState, 'Flrr_volr', begr, endr, input=volr, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport1d(exportState, 'Flrr_volrmch', begr, endr, input=volrmch, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( flds_r2l_stream_channel_depths ) then
       call state_setexport1d(exportState, 'Sr_tdepth', begr, endr, input=tdepth, do_area_correction=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call state_setexport1d(exportState, 'Sr_tdepth_max', begr, endr, input=tdepth_max, do_area_correction=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine export_fields

  !===============================================================================
  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    integer,          optional, intent(in)    :: ungridded_lbound
    integer,          optional, intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname='(rof_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================
  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: dbrc
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(rof_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                     gridToFieldMap=(/2/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
             end if
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(rof_import_export:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================
  subroutine state_getimport1d(state, fldname, begr, endr, area, output, do_area_correction, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr
    integer             , intent(in)    :: endr
    real(r8)            , intent(in)    :: area(begr:endr)
    logical             , intent(in)    :: do_area_correction
    real(r8)            , intent(out)   :: output(begr:endr)
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)            :: lfield
    integer                     :: nr
    real(R8), pointer           :: fldptr(:)
    character(len=*), parameter :: subname='(rof_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! get field pointer
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine output array and scale by unit convertsion
    if (do_area_correction) then
       fldptr(:) = fldptr(:) * med2mod_areacor(:)
    end if
    do nr = begr,endr
       output(nr) = fldptr(nr-begr+1) * area(nr)*0.001_r8
    end do

    ! check for nans
    call check_for_nans(fldptr, trim(fldname), 'import')

  end subroutine state_getimport1d

  !===============================================================================
  subroutine state_getimport2d(state, fldname, begr, endr, ntracers, &
     area, output, do_area_correction, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr
    integer             , intent(in)    :: endr
    integer             , intent(in)    :: ntracers
    real(r8)            , intent(in)    :: area(begr:endr)
    logical             , intent(in)    :: do_area_correction
    real(r8)            , intent(out)   :: output(begr:endr,ntracers)
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)            :: lfield
    integer                     :: nr, nt
    real(r8), pointer           :: fldptr2d(:,:)
    character(len=*), parameter :: subname='(rof_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! get field pointer
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine output array and scale by unit convertsion
    if (do_area_correction) then
      do nt = 1,ntracers
         fldptr2d(nt,:) = fldptr2d(nt,:) * med2mod_areacor(:)
      end do
    end if
    do nt = 1,ntracers
       do nr = begr,endr
          output(nr,nt) = fldptr2d(nt,nr-begr+1) * area(nr)*0.001_r8
       end do
    end do

    ! check for nans
    call check_for_nans(fldptr2d, trim(fldname), 'import')

  end subroutine state_getimport2d

  !===============================================================================
  subroutine state_setexport1d(state, fldname, begr, endr, input, do_area_correction, rc)

    ! ----------------------------------------------
    ! Map input array to export state field
    ! ----------------------------------------------

    use shr_const_mod, only : fillvalue=>shr_const_spval

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr
    integer             , intent(in)    :: endr
    real(r8)            , intent(in)    :: input(begr:endr)
    logical             , intent(in)    :: do_area_correction
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)            :: lfield
    integer                     :: nr
    real(R8), pointer           :: fldptr(:)
    character(len=*), parameter :: subname='(rof_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! get field pointer
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set fldptr values to input array
    fldptr(:) = 0._r8
    do nr = begr,endr
       fldptr(nr-begr+1) = input(nr)
    end do
    if (do_area_correction) then
       fldptr(:) = fldptr(:) * mod2med_areacor(:)
    end if

    ! check for nans
    call check_for_nans(fldptr, trim(fldname), 'export')

  end subroutine state_setexport1d

  !===============================================================================
  subroutine state_setexport2d(state, fldname, begr, endr, ntracers, input, do_area_correction, rc)

    ! ----------------------------------------------
    ! Map input array to export state field
    ! ----------------------------------------------

    use shr_const_mod, only : fillvalue=>shr_const_spval

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr
    integer             , intent(in)    :: endr
    integer             , intent(in)    :: ntracers
    real(r8)            , intent(in)    :: input(begr:endr,ntracers)
    logical             , intent(in)    :: do_area_correction
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)            :: lfield
    integer                     :: nr, nt
    real(R8), pointer           :: fldptr2d(:,:)
    character(len=*), parameter :: subname='(rof_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! get field pointer
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set fldptr values to input array
    fldptr2d(:,:) = 0._r8
    do nt = 1,ntracers
       do nr = begr,endr
          fldptr2d(nt,nr-begr+1) = input(nr,nt)
       end do
       if (do_area_correction) then
          fldptr2d(nt,:) = fldptr2d(nt,:) * mod2med_areacor(:)
       end if
    end do

    ! check for nans
    call check_for_nans(fldptr2d, trim(fldname), 'export')

  end subroutine state_setexport2d

  !===============================================================================

  subroutine check_for_nans_1d(array, fname, type)
     ! Check if any input from mediator or output to mediator is NaN

     ! input/output variables
     real(r8)         , pointer    :: array(:)
     character(len=*) , intent(in) :: fname
     character(len=*) , intent(in) :: type
     !-------------------------------------------------------------------------------

     if (any(isnan(array))) then
        write(iulog,*) "NaN found in field ", trim(fname)
        write(iulog,*) '# of NaNs = ', count(isnan(array))
        write(iulog,*) 'Which are NaNs = ', isnan(array)
        call shr_sys_abort(' ERROR: One or more of the '//trim(type)//' values are NaN ' )
     end if
  end subroutine check_for_nans_1d

  subroutine check_for_nans_2d(array, fname, type)
     ! Check if any input from mediator or output to mediator is NaN

     ! input/output variables
     real(r8)         , pointer    :: array(:,:)
     character(len=*) , intent(in) :: fname
     character(len=*) , intent(in) :: type
     !-------------------------------------------------------------------------------

     if (any(isnan(array))) then
        write(iulog,*) "NaN found in field ", trim(fname)
        write(iulog,*) '# of NaNs = ', count(isnan(array))
        write(iulog,*) 'Which are NaNs = ', isnan(array)
        call shr_sys_abort(' ERROR: One or more of the '//trim(type)//' values are NaN ' )
     end if
  end subroutine check_for_nans_2d

  !===============================================================================
  logical function fldchk(state, fldname)
    ! ----------------------------------------------
    ! Determine if field with fldname is in the input state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag)   :: itemFlag
    ! ----------------------------------------------
    call ESMF_StateGet(state, trim(fldname), itemFlag)
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       fldchk = .true.
    else
       fldchk = .false.
    endif
  end function fldchk

end module rof_import_export
