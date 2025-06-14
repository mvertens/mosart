module mosart_physics

   !-----------------------------------------------------------------------
   ! Description: core code of MOSART.
   ! Contains routines for solving diffusion wave and update the state of
   ! hillslope, subnetwork and main channel variables
   ! Developed by Hongyi Li, 12/29/2011.
   !-----------------------------------------------------------------------

   use shr_kind_mod      , only : r8 => shr_kind_r8
   use shr_const_mod     , only : SHR_CONST_REARTH, SHR_CONST_PI
   use shr_sys_mod       , only : shr_sys_abort
   use mosart_vars       , only : iulog, barrier_timers, mpicom_rof, bypass_routing_option,        &
                                  debug_mosart
   use mosart_data       , only : Tctl, TUnit, TRunoff, TPara, ctl
   use perf_mod          , only : t_startf, t_stopf
   use nuopc_shr_methods , only : chkerr
   use ESMF              , only : ESMF_FieldGet, ESMF_FieldSMM, ESMF_Finalize, &
                                  ESMF_SUCCESS, ESMF_END_ABORT, ESMF_TERMORDER_SRCSEQ

   implicit none
   private

   public  :: Euler
   public  :: mosart_physics_restart

   private :: hillsloperouting
   private :: updatestate_hillslope

   private :: subnetworkrouting
   private :: updatestate_subnetwork

   private :: mainchannelrouting
   private :: updatestate_mainchannel

   private :: CRVRMAN ! Function for calculating channel velocity according to Manning's equation.
   private :: CREHT   ! Function for overland from hillslope into the sub-network channels
   private :: GRMR    ! Function for estimate wetted channel area
   private :: GRHT    ! Function for estimating water depth assuming rectangular channel
   private :: GRPT    ! Function for estimating wetted perimeter assuming rectangular channel
   private :: GRRR    ! Function for estimating hydraulic radius
   private :: GRPR    ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain

   real(r8), parameter :: TINYVALUE = 1.0e-14_r8    ! double precision variable has a significance of about 16 decimal digits
   real(r8), parameter :: limiter_advec = 0.9999_r8 ! limit advection to maximum close to 1 of available storage(+inflow),
                                                    ! courant number cn=1 leads to numerical precision issues (small negative values)
   real(r8), parameter :: SLOPE1def = 0.1_r8     ! here give it a small value in order to avoid the abrupt change of hydraulic radidus etc.
   real(r8)            :: sinatanSLOPE1defr      ! 1.0/sin(atan(slope1))

   character(*), parameter :: u_FILE_u = &
        __FILE__

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

   subroutine Euler(rc)

      ! solve the ODEs with Euler algorithm

      ! Arguments
      integer, intent(out) :: rc

      ! Local variables
      integer           :: nt, nr, nt_nonh2o, m, k, unitUp, cnt, ier   !local index
      real(r8)          :: localDeltaT
      real(r8)          :: negchan
      real(r8), pointer :: src_eroutUp(:,:)
      real(r8), pointer :: dst_eroutUp(:,:)
      real(r8),allocatable :: temp_erout(:)

      rc = ESMF_SUCCESS

      associate( &
           ! ctl
           begr        => ctl%begr,            & ! local start index
           endr        => ctl%endr,            & ! local end index
           ntracers    => ctl%ntracers_tot,    & ! total number of tracers
           ntracers_nonH2O=>ctl%ntracers_nonh2o,& ! number of non-standard H2O tracers
           nt_liq      => ctl%nt_liq,          & ! index for liquid water
           nt_ice      => ctl%nt_ice,          & ! index for ice

           ! tctl
           deltat      => tctl%DeltaT,         & ! Time step [s]
           DLevelH2R   => tctl%DLevelH2R,      & ! The base number of channel routing sub-time-steps within one hillslope routing step.

           ! Tunit
           euler_calc  => Tunit%euler_calc,    & ! flag to calculate tracers in euler
           area        => Tunit%area,          & ! area of local cell [m2]
           frac        => Tunit%frac,          & ! fraction of cell included in the study area, [-]
           mask        => Tunit%mask,          & ! mask for cell 0=null, 1=land with dnID, 2=outlet
           numDT_r     => Tunit%numDT_r,       & ! for a main reach, the number of sub-time-steps needed for numerical stability
           numDT_t     => Tunit%numDT_t,       & ! for a subnetwork reach, the number of sub-time-steps needed for numerical stability

           ! hillsloope
           !! states
           wh          => TRunoff%wh,          & ! storage of surface water, [m] or tracer [kg/m2, mol/m2, or similar]
           dwh         => TRunoff%dwh,         & ! change of water storage, [m/s] or tracer [kg/m2/s, mol/m2/s, or similar]
           yh          => TRunoff%yh,          & ! depth of surface water, [m]
           wsat        => Trunoff%wsat,        & ! storage of surface water within saturated area at hillslope [m]    - NOT USED
           wunsat      => Trunoff%wunsat,      & ! storage of surface water within unsaturated area at hillslope [m]  - NOT USED
           qhorton     => Trunoff%qhorton,     & ! Infiltration excess runoff generated from hillslope, [m/s]         - NOT_USED
           qdunne      => Trunoff%qdunne,      & ! Saturation excess runoff generated from hillslope, [m/s]           - NOT_USED
           qsur        => Trunoff%qsur,        & ! Surface runoff generated from hillslope, [m/s]
           qsub        => Trunoff%qsub,        & ! Subsurface runoff generated from hillslope, [m/s]
           qgwl        => Trunoff%qgwl,        & ! gwl runoff term from glacier, wetlands and lakes, [m/s]
           !! fluxes
           ehout       => Trunoff%ehout,       & ! overland flow from hillslope into the sub-channel, water [m/s] or tarcers [kg/m2/s, mol/m2/s, or similar]

           ! subnetwork channel
           !! states
           tarea       => Trunoff%tarea,       & ! area of channel water surface, [m2]
           wt          => Trunoff%wt,          & ! storage of surface water, [m3] or tracer [kg, mol, or similar]
           dwt         => Trunoff%dwt,         & ! change of water storage, [m3/s] or tracer [kg/s, mol/s, or similar]
           yt          => Trunoff%yt,          & ! water depth, [m]
           mt          => Trunoff%mt,          & ! cross section area, [m2]
           pt          => Trunoff%pt,          & ! wetness perimeter, [m]
           vt          => Trunoff%vt,          & ! flow velocity, [m/s]
           rt          => TRunoff%rt,          & ! hydraulic radii, [m]
           !! fluxes
           etin        => Trunoff%etin,        & ! lateral water [m3/s] or tracer [kg/s, mol/s or similar] inflow from hillslope
           etout       => Trunoff%etout,       & ! water [ m3/s] or tracer [kg/s, mol/s or similar] discharge from sub-network into the main reach

           ! main channel
           !! states
           rarea       => Trunoff%rarea,       & ! area of channel water surface, [m2]
           wr          => Trunoff%wr,          & ! storage of surface water, [m3] or tracer [kg, mol, or similar]
           dwr         => Trunoff%dwr,         & ! change of water storage, [m3] or tracer [kg/s, mol/s or similar]
           yr          => Trunoff%yr,          & ! water depth. [m]
           mr          => Trunoff%mr,          & ! cross section area, [m2]
           rr          => Trunoff%rr,          & ! hydraulic radius, [m]
           pr          => Trunoff%pr,          & ! wetness perimeter, [m]
           vr          => Trunoff%vr,          & ! flow velocity, [m/s]
           !! exchange fluxes
           erlateral   => Trunoff%erlateral,   & ! lateral water flow [m3/s] or tracer [kg/s, mol/s or similar] from hillslope
           erin        => Trunoff%erin,        & ! inflow from upstream links, water [m3/s] or tracer [kg/s, mol/s or similar]
           erout       => Trunoff%erout,       & ! outflow into downstream links, water [m3/s] or tracer [kg/s, mol/s or similar]
           erout_prev  => Trunoff%erout_prev,  & ! outflow into downstream links from previous timestep, water [m3/s] or tracer [kg/s, mol/s or similar]
           eroutUp     => Trunoff%eroutUp,     & ! outflow sum of upstream gridcells, instantaneous water [m3/s] or tracer [kg/s, mol/s or similar]
           eroutUp_avg => Trunoff%eroutUp_avg, & ! outflow sum of upstream gridcells, average water [m3/s] or tracer [kg/s, mol/s or similar]
           erlat_avg   => Trunoff%erlat_avg,   & ! erlateral average water [m3/s] or tracer [kg/s, mol/s or similar]
           flow        => Trunoff%flow         & ! streamflow from the outlet of the reach, water [m3/s] or tracer [kg/s, mol/s or similar]
      )

      allocate(temp_erout(ntracers))

      !------------------
      ! hillslope
      !------------------

      call t_startf('mosartr_hillslope')
      do nt = 1,ntracers ! could be simplified to nt=1 and no if(euler_calc)
         if (euler_calc(nt)) then
            do nr = begr,endr
               if (mask(nr) > 0) then
                  ! Liquid water flow - calculate time differences
                  call hillslopeRouting(nr, DeltaT, yh(nr,nt), wh(nr,nt), qsur(nr,nt),             & ! input
                                        ehout(nr,nt), dwh(nr,nt))                                    ! ouput

                  ! Advect tracers with water flow - calculate time differences
                  if (ntracers_nonH2O > 0 .and. nt == nt_liq) then
                    do nt_nonh2o = nt_ice+1,nt_ice+ntracers_nonH2O
                      call HSR_advect_tracers(nr,DeltaT,wh(nr,nt_liq),ehout(nr,nt_liq),            & ! water mass advection info
                                              wh(nr,nt_nonh2o),qsur(nr,nt_nonh2o),                 & ! nonH2O storage and input
                                              ehout(nr,nt_nonh2o),dwh(nr,nt_nonh2o))                 ! output (mass outflow & mass storage change)
                      wh(nr,nt_nonh2o)  = wh(nr,nt_nonh2o) + dwh(nr,nt_nonh2o)*DeltaT                ! update mass storage [here kg/m2 or mol/m2]
                      etin(nr,nt_nonh2o)= (-ehout(nr,nt_nonh2o)+qsub(nr,nt_nonh2o)) * area(nr) * frac(nr)  ! update lateral + surf + subsurf mass inflow [kg/s or mol/s]
                      if (debug_mosart >= 2) then
                        if (wh(nr,nt_nonh2o) <0._r8)then
                            write(iulog,*) 'DEBUG: HSR Warning: Negative storage found! ',         &
                                     nr,wh(nr,nt_nonh2o),dwh(nr,nt_nonh2o)*DeltaT,qsur(nr,nt_nonh2o)
                        endif
                      endif
                    enddo
                  endif ! end nonH2O

                  ! Updating water storage and tracers to next time step
                  wh(nr,nt) = wh(nr,nt) + dwh(nr,nt) * DeltaT
                  call UpdateState_hillslope(wh(nr,nt), yh(nr,nt))
                  etin(nr,nt) = (-ehout(nr,nt) + qsub(nr,nt)) * area(nr) * frac(nr)
                  if (debug_mosart >= 1) then
                    ! check for negative hillslope storage
                    if(wh(nr,1) < 0._r8) then
                      write(iulog,*) 'DEBUG: Negative hillslope storage! ', nr, wh(nr,1)
                      !call shr_sys_abort('mosart: negative hillslope storage')
                    end if
                  endif
               endif
            end do
         endif
      end do
      call t_stopf('mosartr_hillslope')

      !------------------
      ! subnetwork and channel
      !------------------

      call ESMF_FieldGet(Tunit%srcfield, farrayPtr=src_eroutUp, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldGet(Tunit%dstfield, farrayPtr=dst_eroutUp, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      src_eroutUp(:,:) = 0._r8
      dst_eroutUp(:,:) = 0._r8

      Trunoff%flow(:,:) = 0._r8
      Trunoff%erout_prev(:,:) = 0._r8
      Trunoff%eroutup_avg(:,:) = 0._r8
      Trunoff%erlat_avg(:,:) = 0._r8
      negchan = 9999.0_r8

      do m = 1,DLevelH2R

         ! accumulate/average erout at prior timestep (used in eroutUp calc) for budget analysis
         do nt=1,ntracers
           do nr=begr,endr
              Trunoff%erout_prev(nr,nt) = Trunoff%erout_prev(nr,nt) + erout(nr,nt)
           end do
         end do

         !------------------
         ! subnetwork
         !------------------

         call t_startf('mosartr_subnetwork')
         erlateral(:,:) = 0._r8
         do nt=1,ntracers
            if (euler_calc(nt)) then ! could be simplified to nt=1
               do nr = begr,endr
                  if (mask(nr) > 0) then
                     localDeltaT = DeltaT/DLevelH2R/numDT_t(nr)
                     do k = 1,numDT_t(nr)
                        ! Liquid water flow - calculate time differences
                        call subnetworkRouting(nr, localDeltaT, rt(nr,nt), mt(nr,nt),              & ! input
                                               wt(nr,nt), etin(nr,nt),                             & ! input
                                               etout(nr,nt), vt(nr,nt), dwt(nr,nt))                  ! output

                        ! Advect tracers with water flow - calculate time differences
                        if (ntracers_nonH2O > 0 .and. nt == nt_liq) then
                          do nt_nonh2o = nt_ice+1,nt_ice+ntracers_nonH2O
                            call SNR_advect_tracers(nr,localDeltaT,wt(nr,nt_liq),etout(nr,nt_liq), & ! water mass advection info
                                                    wt(nr,nt_nonh2o),etin(nr,nt_nonh2o),           & ! nonH2O storage and input
                                                    etout(nr,nt_nonh2o),dwt(nr,nt_nonh2o))           ! output (mass outflow & mass storage change)
                            wt(nr,nt_nonh2o)  = wt(nr,nt_nonh2o) + dwt(nr,nt_nonh2o)*localDeltaT     ! update mass storage [here: kg or mol]
                            erlateral(nr,nt_nonh2o)= erlateral(nr,nt_nonh2o) - etout(nr,nt_nonh2o)   ! update lateral + surf + subsurf mass inflow [kg/s or mol/s]
                            if (debug_mosart >= 2) then
                              if (wt(nr,nt_nonh2o) < 0._r8)then
                               write(iulog,*) 'DEBUG: SNR Warning: Negative storage found! ',      &
                                          nr,wt(nr,nt_nonh2o),dwt(nr,nt_nonh2o)*localDeltaT,       &
                                          etin(nr,nt_nonh2o),etout(nr,nt_nonh2o)
                              endif
                            endif
                          enddo
                        endif ! end nonH2O

                        ! Updating water storage and tracers to next time step
                        wt(nr,nt) = wt(nr,nt) + dwt(nr,nt) * localDeltaT
                        if (debug_mosart >= 1) then
                          ! check for negative subnetwork storage
                          if(wt(nr,1) < 0._r8) then
                            write(iulog,*) 'DEBUG: Negative subnetwork storage! ', nr, wt(nr,1)
                            !call shr_sys_abort('mosart: negative subnetwork storage')
                          end if
                        endif
                        call UpdateState_subnetwork(nr, wt(nr,nt),                                 &     ! input
                             mt(nr,nt), yt(nr,nt), pt(nr,nt), rt(nr,nt))                                 ! output
                        erlateral(nr,nt) = erlateral(nr,nt) - etout(nr,nt)
                     end do ! numDT_t
                     erlateral(nr,nt) = erlateral(nr,nt) / numDT_t(nr)
                     if (ntracers_nonH2O > 0 .and. nt == nt_liq) then
                        do nt_nonh2o = nt_ice+1,nt_ice+ntracers_nonH2O
                          erlateral(nr,nt_nonh2o) = erlateral(nr,nt_nonh2o) / numDT_t(nr)
                        end do
                     endif
                  endif
               end do ! nr
            endif  ! euler_calc
         end do ! nt
         call t_stopf('mosartr_subnetwork')

         !------------------
         ! upstream interactions
         !------------------

         if (barrier_timers) then
            call t_startf('mosartr_SMeroutUp_barrier')
            call mpi_barrier(mpicom_rof,ier)
            call t_stopf('mosartr_SMeroutUp_barrier')
         endif

         call t_startf('mosartr_SMeroutUp')

         !--- copy erout into src_eroutUp ---
         eroutUp = 0._r8
         src_eroutUp(:,:) = 0._r8
         cnt = 0
         do nr = begr,endr
            cnt = cnt + 1
            do nt = 1,ntracers
               src_eroutUp(nt,cnt) = erout(nr,nt)
            enddo
         enddo

         ! --- map src_eroutUp to dst_eroutUp
         call ESMF_FieldSMM(TUnit%srcfield, TUnit%dstField, TUnit%rh_eroutUp, termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         !--- copy mapped eroutUp to TRunoff ---
         cnt = 0
         do nr = begr,endr
            cnt = cnt + 1
            do nt = 1,ntracers
               eroutUp(nr,nt) = dst_eroutUp(nt,cnt)
            enddo
         enddo

         call t_stopf('mosartr_SMeroutUp')

         Trunoff%eroutup_avg = Trunoff%eroutup_avg + Trunoff%eroutUp
         Trunoff%erlat_avg   = Trunoff%erlat_avg + Trunoff%erlateral

         !------------------
         ! main channel routing
         !------------------

         call t_startf('mosartr_chanroute')

         do nt = 1,ntracers
            if (euler_calc(nt)) then
               do nr = begr,endr
                  if(mask(nr) > 0) then
                     localDeltaT = DeltaT/DLevelH2R/numDT_r(nr)
                     temp_erout = 0._r8
                     do k = 1,numDT_r(nr)
                        ! Liquid water flow - calculate time differences
                        call mainchannelRouting(nr, localDeltaT,                                   & ! input
                             eroutUp(nr,nt), erlateral(nr,nt),                                     & ! input
                             wr(nr,nt),mr(nr,nt), rr(nr,nt), qgwl(nr,nt),                          & ! input
                             erin(nr,nt), erout(nr,nt), vr(nr,nt), dwr(nr,nt))                       ! output

                        ! Advect tracers with water flow - calculate time differences
                        if (ntracers_nonH2O > 0 .and. nt == nt_liq) then
                          do nt_nonh2o = nt_ice+1,nt_ice+ntracers_nonH2O
                            call MCR_advect_tracers(nr, localDeltaT,eroutUp(nr,nt_nonh2o),         & ! input
                                                    erlateral(nr,nt_nonh2o),                       & ! input
                                                    wr(nr,nt_liq), mr(nr,nt_liq), vr(nr,nt_liq),   & ! input
                                                    rr(nr,nt_liq),                                 & ! input
                                                    wr(nr,nt_nonh2o), qgwl(nr,nt_nonh2o),          & ! input
                                                    erin(nr,nt_nonh2o),erout(nr,nt_nonh2o),        & ! output
                                                    dwr(nr,nt_nonh2o))                               ! output
                            wr(nr,nt_nonh2o)  = wr(nr,nt_nonh2o) + dwr(nr,nt_nonh2o)*localDeltaT     ! update mass storage [here: kg or mol]
                            temp_erout(nt_nonh2o) = temp_erout(nt_nonh2o) + erout(nr,nt_nonh2o)
                            if (debug_mosart >= 2) then
                              if (wr(nr,nt_nonh2o) <0._r8)then
                                write(iulog,*) 'DEBUG: MCR Warning: Negative storage found! ',     &
                                        nr,wr(nr,nt_nonh2o),dwr(nr,nt_nonh2o)*localDeltaT
                              endif
                            endif
                          enddo
                        endif ! end nonH2O

                        ! Updating water storage and tracers to next time step
                        wr(nr,nt) = wr(nr,nt) + dwr(nr,nt) * localDeltaT
                        if (debug_mosart >= 1) then
                          ! check for negative channel storage
                          if(wr(nr,1) < 0._r8) then
                            write(iulog,*) 'DEBUG: Negative channel storage! ', nr, wr(nr,1)
                            !call shr_sys_abort('mosart: negative channel storage')
                          end if
                        end if

                        call UpdateState_mainchannel(nr, wr(nr,nt), &    ! input
                             mr(nr,nt), yr(nr,nt), pr(nr,nt), rr(nr,nt)) ! output

                        ! erout here might be inflow to some downstream subbasin,
                        ! so treat it differently than erlateral
                        temp_erout(nt) = temp_erout(nt) + erout(nr,nt)
                     end do
                     temp_erout(nt) = temp_erout(nt) / numDT_r(nr)
                     erout(nr,nt) = temp_erout(nt)
                     Trunoff%flow(nr,nt) = Trunoff%flow(nr,nt) - erout(nr,nt)
                     if (ntracers_nonH2O > 0 .and. nt == nt_liq) then
                       do nt_nonh2o = nt_ice+1,nt_ice+ntracers_nonH2O
                         erout(nr,nt_nonh2o) = temp_erout(nt_nonh2o) /numDT_r(nr)
                         Trunoff%flow(nr,nt_nonh2o) = Trunoff%flow(nr,nt_nonh2o) - erout(nr,nt_nonh2o)
                       enddo
                     endif
                  endif
               end do ! nr
            endif  ! euler_calc
         end do ! nt
         negchan = min(negchan, minval(wr(:,:)))

         call t_stopf('mosartr_chanroute')
      end do ! DLevelH2R

      if (debug_mosart >= 1) then
        ! check for negative channel storage
        if (negchan < 0._r8) then
           write(iulog,*) 'DEBUG: Warning: Negative channel storage found! ',negchan
           ! call shr_sys_abort('mosart: negative channel storage')
        endif
      endif
      flow        = flow        / DLevelH2R
      erout_prev  = erout_prev  / DLevelH2R
      eroutup_avg = eroutup_avg / DLevelH2R
      erlat_avg   = erlat_avg   / DLevelH2R

    end associate

   end subroutine Euler

   !-----------------------------------------------------------------------

   subroutine hillslopeRouting(nr, DeltaT, yh, wh, qsur, ehout, dwh)
      !  Hillslope routing considering uniform runoff generation across hillslope

      ! Arguments
      integer , intent(in)  :: nr     ! index number to use
      real(r8), intent(in)  :: DeltaT ! Time step in seconds
      real(r8), intent(in)  :: yh     ! depth of surface water, [m]
      real(r8), intent(in)  :: wh     ! storage of surface water, [m]
      real(r8), intent(in)  :: qsur   ! Surface runoff generated from hillslope, [m/s]
      real(r8), intent(out) :: ehout  ! overland flow from hillslope into the sub-channel, [m/s]
      real(r8), intent(out) :: dwh    ! change of water storage, [m/s]

      ehout = -CREHT(TUnit%hslpsqrt(nr), TUnit%nh(nr), TUnit%Gxr(nr), yh)
      if(ehout < 0._r8 .and. (wh + (qsur + ehout) * DeltaT) < TINYVALUE) then
         ! limit the amount of outflow to limiter_advec fraction and prevent from going negative
         ehout = -limiter_advec*(qsur + max(0._r8,wh-epsilon(1._r8)) / DeltaT)
      end if
      dwh = (qsur + ehout)

   end subroutine hillslopeRouting

   !-----------------------------------------------------------------------

   subroutine HSR_advect_tracers(nr,DeltaT,wh_liq,ehout_liq,wh,qsur,ehout,dwh)
     ! Hillslope routing advection of tracers

     ! Arguments
     integer, intent(in) :: nr        ! index number to use
     real(r8),intent(in) :: DeltaT    ! Time step in seconds
     real(r8),intent(in) :: wh_liq    ! storage of surface water [m]
     real(r8),intent(in) :: ehout_liq ! overland flow from hillslope into the sub-channel [m/s]
     real(r8),intent(in) :: wh        ! mass storage [kg/m2 or mol/m2]
     real(r8),intent(in) :: qsur      ! surface mass inflow [kg/m2/s or mol/m2/s]
     real(r8),intent(out):: ehout     ! overland mass flow from hillslopes into sub-channels [kg/m2/s or mol/m2/s]
     real(r8),intent(out):: dwh       ! change of mass storage [kg/m2/s or mol/m2/s]

     if (abs(-1._r8/limiter_advec*ehout_liq*DeltaT - wh_liq) < TINYVALUE) then
        ehout= -limiter_advec*(qsur + max(0._r8,wh-epsilon(1.))/DeltaT)
     else
        ehout = ehout_liq/(wh_liq+epsilon(1._r8)) * wh
     endif
     dwh   = qsur + ehout

   end subroutine HSR_advect_tracers

   !-----------------------------------------------------------------------

   subroutine subnetworkRouting(nr, DeltaT, rt, mt, wt, etin, etout, vt, dwt)
      !  subnetwork channel routing

      ! Arguments
      integer , intent(in)  :: nr    ! index number to use
      real(r8), intent(in)  :: DeltaT! Time step [s]
      real(r8), intent(in)  :: rt    ! hydraulic radius, [m]
      real(r8), intent(in)  :: mt    ! cross section area, [m2]
      real(r8), intent(in)  :: wt    ! storage of surface water, [m3]
      real(r8), intent(in)  :: etin  ! lateral inflow from hillslope
      real(r8), intent(out) :: etout ! discharge from sub-network into the main reach, [m3/s]
      real(r8), intent(out) :: vt    ! flow velocity, [m/s]
      real(r8), intent(out) :: dwt   ! change of water storage, [m3]

      if(TUnit%tlen(nr) <= TUnit%hlen(nr)) then ! if no tributaries, not subnetwork channel routing
         etout = -etin
      else
         vt = CRVRMAN(TUnit%tslpsqrt(nr), TUnit%nt(nr), rt)
         etout = -vt * mt
         if (wt + (etin + etout) * DeltaT < TINYVALUE) then
            etout = -limiter_advec*(etin + max(0._r8,wt-epsilon(1._r8))/DeltaT)
            if (mt > 0._r8) then
               vt = -etout/mt
            end if
         end if
      end if
      dwt = etin + etout

      if (debug_mosart >= 3) then
        ! check stability
        if (vt < 0._r8 .or. vt > 30._r8) then
           write(iulog,*) "DEBUG: Numerical error in subnetworkRouting flow velocity, ", nr,vt
        end if
      endif
   end subroutine subnetworkRouting

   !-----------------------------------------------------------------------

   subroutine SNR_advect_tracers(nr,DeltaT, wt_liq,etout_liq,wt,etin,etout,dwt)
     ! Subnetwork channel routing tracer advection

     ! Arguments
     integer, intent(in) :: nr        ! index number to use
     real(r8),intent(in) :: DeltaT    ! Time step [s]
     real(r8),intent(in) :: wt_liq    ! storage of surface water [m3]
     real(r8),intent(in) :: etout_liq ! sub-channel flow [m3/s]
     real(r8),intent(in) :: wt        ! mass storage [kg or mol]
     real(r8),intent(in) :: etin      ! mass inflow [kg/s or mol/s]
     real(r8),intent(out):: etout     ! mass flow from sub-channels into main-reach [kg/s or mol/s]
     real(r8),intent(out):: dwt       ! change of mass storage [kg/s or mol/s]

     if (TUnit%tlen(nr) <= TUnit%hlen(nr)) then ! if no tributaries, no subnetwork channel routing
       etout = -etin
     else
        if (-(1._r8/limiter_advec*etout_liq*DeltaT) >= wt_liq) then ! not completely coherent with subchannel routing
           etout = -limiter_advec*(etin + max(0._r8,wt-epsilon(1._r8))/DeltaT)
        else
           etout = etout_liq/(wt_liq+epsilon(1._r8)) * wt
        endif
     endif
     dwt = etin + etout
   end subroutine SNR_advect_tracers

   !-----------------------------------------------------------------------

   subroutine mainchannelRouting(nr, DeltaT, eroutUp, erlateral, wr, mr, rr, qgwl, &
        erin, erout, vr, dwr)

      !  classic kinematic wave routing method

      ! Arguments
      integer , intent(in)  :: nr        ! index number to use
      real(r8), intent(in)  :: DeltaT    ! time step [s]
      real(r8), intent(in)  :: eroutUp   ! outflow sum of upstream gridcells, instantaneous (m3/s)
      real(r8), intent(in)  :: erlateral ! lateral flow from hillslope [m3/s]
      real(r8), intent(in)  :: wr        ! storage of surface water, [m3]
      real(r8), intent(in)  :: mr        ! cross section area, [m2]
      real(r8), intent(in)  :: rr        ! hydraulic radius, [m]
      real(r8), intent(in)  :: qgwl      ! gwl runoff term from glacier, wetlands and lakes, [m/s]
      real(r8), intent(out) :: erin      ! inflow from upstream links, [m3/s]
      real(r8), intent(out) :: erout     ! outflow into downstream links, [m3/s]
      real(r8), intent(out) :: vr        ! flow velocity, [m/s]
      real(r8), intent(out) :: dwr       ! change in water storage, [m3]

      ! Local variables
      real(r8) :: temp_gwl

      associate( &
           roughl     => Tunit%nr  ,       & ! manning's roughness of the main reach
           area       => Tunit%area,       & ! area of the gridcell [m2]
           areaTotal2 => Tunit%areaTotal2, & ! computed total upstream drainage area, [m2]
           frac       => Tunit%frac,       & ! fraction of cell included in the study area, [-]
           rlen       => Tunit%rlen,       & ! length of main river reach, [m]
           rwidth     => Tunit%rwidth,     & ! bankfull width of main reach, [m]
           twidth     => Tunit%twidth,     & ! bankfull width of the sub-reach, [m]
           rslpsqrt   => TUnit%rslpsqrt    & ! sqrt of slope of main river reach, [-]
         )

      ! estimate the inflow from upstream units
      erin = 0._r8
      erin = erin - eroutUp

      ! estimate the outflow and flow velocity
      if (rlen(nr) <= 0._r8) then ! no river network, no channel routing
         vr = 0._r8
         erout = -erin - erlateral
      else
        if(areaTotal2(nr)/(rwidth(nr)*rlen(nr)) > 1e6_r8) then ! if total drainage area is much larger than channel surface area (large river)
            erout = -erin - erlateral
        else
            vr = CRVRMAN(rslpsqrt(nr), roughl(nr), rr)
            erout = -vr * mr
            if (-erout > TINYVALUE .and. wr + (erlateral + erin + erout) * DeltaT < TINYVALUE) then
               erout = -(erlateral + erin + wr / DeltaT)
               if (mr > 0._r8) then ! why should cross section area become <= zero
                  vr = -erout / mr
               end if
            end if
        end if
      end if

      temp_gwl = qgwl * area(nr) * frac(nr)
      dwr = erlateral + erin + erout + temp_gwl

      if ((wr/DeltaT + dwr) < 0._r8 .and. (trim(bypass_routing_option)/='none') ) then
         write(iulog,*) 'DEBUG: mosart: ERROR main channel going negative: ', nr
         write(iulog,*) DeltaT, wr, wr/DeltaT, dwr, temp_gwl
         write(iulog,*) ' '
      endif

      if (debug_mosart >= 1) then
        ! check for stability
        if(vr < 0._r8 .or. vr > 30._r8) then
           write(iulog,*) "DEBUG: Numerical error in Routing_KW flow velocity, ", nr,vr
        end if

        ! check for negative wr
        if(wr > 1._r8 .and. (wr/DeltaT + dwr)/wr < 0._r8) then
           write(iulog,*) 'DEBUG: negative wr!', wr, dwr, temp_gwl, DeltaT
        !       stop
        end if
      endif
      end associate
   end subroutine MainchannelRouting

   !-----------------------------------------------------------------------

   subroutine MCR_advect_tracers(nr,DeltaT, eroutUp, erlateral, wr_liq, mr, vr, rr, wr, qgwl, erin, erout, dwr)
     ! Main channel routing tracer advection
     !Arguments
     integer,  intent(in)  :: nr        ! index number to use
     real(r8), intent(in)  :: DeltaT    ! time step [s]
     real(r8), intent(in)  :: eroutUp   ! mass outflow sum from upstream gridcells, instanatneous [kg/s or mol/s]
     real(r8), intent(in)  :: erlateral ! lateral mass flow from hillslope [kg/s or mol/s]
     real(r8), intent(in)  :: wr_liq    ! water storage [m3]
     real(r8), intent(in)  :: mr        ! cross section area [m2]
     real(r8), intent(in)  :: vr        ! flow velocity [m/s]
     real(r8), intent(in)  :: rr        ! hydraulic radius [m]
     real(r8), intent(in)  :: wr        ! tracer mass storage [kg or mol]
     real(r8), intent(in)  :: qgwl      ! gwl tracer mass flux from glacier, wetlands and lakes [kg/m2/s or mol/m2/s]
     real(r8), intent(out) :: erin      ! mass inflow from upstream links [kg/s or mol/s]
     real(r8), intent(out) :: erout     ! mass inflow into downstream links [kg/s or mol/s]
     real(r8), intent(out) :: dwr       ! change in mass storage

     associate(area       => Tunit%area,       & ! area of the gridcell [m2]
               areaTotal2 => Tunit%areaTotal2, & ! computed total upstream drainage area, [m2]
               frac       => Tunit%frac,       & ! fraction of cell included in the study area, [-]
               rlen       => Tunit%rlen,       & ! length of main river reach, [m]
               rwidth     => Tunit%rwidth,     & ! bankfull width of main reach, [m]
               roughl     => Tunit%nr  ,       & ! manning's roughness of the main reach
               rslpsqrt   => TUnit%rslpsqrt    & ! sqrt of slope of main river reach [-]
              )

     ! Mass inflow from upstream regions
     erin = 0._r8
     erin = erin - eroutUp

     ! estimate the outflow
     if (rlen(nr) <= 0._r8) then ! no river network, no channel routing
       erout = -erin - erlateral
     else
       if(areaTotal2(nr)/(rwidth(nr)*rlen(nr)) > 1e6_r8) then
         erout = -erin - erlateral
       else
         if (abs(vr - CRVRMAN(rslpsqrt(nr), roughl(nr), rr)) < TINYVALUE) then
           erout = -vr*mr/(wr_liq+epsilon(1._r8)) * wr
         else
           erout = -(erlateral + erin + (wr_liq / DeltaT)/(wr_liq+epsilon(1._r8)) * wr)
         endif
       endif
     endif

     dwr = erlateral + erin + erout + qgwl * area(nr) * frac(nr)

     end associate
   end subroutine MCR_advect_tracers
   !-----------------------------------------------------------------------

   subroutine updateState_hillslope(wh, yh)
     !  update the state variables at hillslope

     ! Arguments
     real(r8), intent(in)  :: wh ! storage of surface water, [m]
     real(r8), intent(out) :: yh ! depth of surface water, [m]

     yh = wh !/ TUnit%area(nr) / TUnit%frac(nr)

   end subroutine updateState_hillslope

   !-----------------------------------------------------------------------

   subroutine updateState_subnetwork(nr, wt, mt, yt, pt, rt)
      !  update the state variables in subnetwork channel

      ! Arguments
      integer , intent(in)  :: nr ! index to use
      real(r8), intent(in)  :: wt ! storage of surface water, [m3]
      real(r8), intent(out) :: mt ! cross section area, [m2]
      real(r8), intent(out) :: yt ! water depth, [m]
      real(r8), intent(out) :: pt ! wetness perimeter, [m]
      real(r8), intent(out) :: rt ! hydraulic radii, [m]

      associate( &
           tlen   => TUnit%tlen,  &
           twidth => Tunit%twidth &
      )

      if(tlen(nr) > 0._r8 .and. wt > 0._r8) then
         mt = GRMR(wt, tlen(nr))
         yt = GRHT(mt, twidth(nr))
         pt = GRPT(yt, twidth(nr))
         rt = GRRR(mt, pt)
      else
         mt = 0._r8
         yt = 0._r8
         pt = 0._r8
         rt = 0._r8
      end if

      end associate
   end subroutine updateState_subnetwork

   !-----------------------------------------------------------------------

   subroutine updateState_mainchannel(nr, wr, mr, yr, pr, rr)
      !  update the state variables in main channel

      ! Arguments
      integer, intent(in)   :: nr  ! index value to use
      real(r8), intent(in)  :: wr  ! storage of surface water, [m3]
      real(r8), intent(out) :: mr  ! cross section area, [m2]
      real(r8), intent(out) :: yr  ! water depth. [m]
      real(r8), intent(out) :: pr  ! wetness perimeter, [m]
      real(r8), intent(out) :: rr  ! hydraulic radius, [m]

      associate( &
           rlen    => Tunit%rlen,    & ! length of main river reach, [m]
           rwidth  => Tunit%rwidth,  & ! bankfull width of main reach, [m]
           rwidth0 => Tunit%rwidth0, & ! total width of the flood plain, [m]
           rdepth  => Tunit%rdepth   & ! bankfull depth of river cross section, [m]
      )

      if(TUnit%rlen(nr) > 0._r8 .and. wr > 0._r8) then
         mr = GRMR(wr, rlen(nr))
         yr = GRHR(mr, rwidth(nr), rwidth0(nr), rdepth(nr))
         pr = GRPR(yr, rwidth(nr), rwidth0(nr), rdepth(nr))
         rr = GRRR(mr, pr)
      else
         mr = 0._r8
         yr = 0._r8
         pr = 0._r8
         rr = 0._r8
      end if

      end associate
   end subroutine updateState_mainchannel

   !-----------------------------------------------------------------------

   function CRVRMAN(sqrtslp_, n_, rr_) result(v_)
      ! Function for calculating channel velocity according to Manning's equation.

      ! Arguments
      real(r8), intent(in) :: sqrtslp_ ! sqrt of average slope of tributaries, [-]
      real(r8), intent(in) :: n_       ! manning's roughness coeff. [s/m^(1/3)]
      real(r8), intent(in) :: rr_      ! hydraulic radius [m]
      real(r8)             :: v_       ! v_ is discharge velocity [m/s]

      if (rr_ <= 0._r8) then
         v_ = 0._r8
      else
         v_ = ((rr_*rr_)**(1._r8/3._r8)) * sqrtslp_ / n_
      end if

   end function CRVRMAN

   !-----------------------------------------------------------------------

   function CREHT(sqrthslp_, nh_, Gxr_, yh_) result(eht_)
      ! Function for overland from hillslope into the sub-network channels

      ! Arguments
      real(r8), intent(in) :: sqrthslp_ ! sqrt of slope of hillslope, [-]
      real(r8), intent(in) :: nh_       ! manning's roughness coeff. [s/m^(1/3)]
      real(r8), intent(in) :: Gxr_      ! drainage density within the cell, [1/m]
      real(r8), intent(in) :: yh_       ! depth of surface water, [m]
      real(r8)             :: eht_      ! velocity, specific discharge [m/s]

      ! Local variables
      real(r8) :: vh_

      vh_ = CRVRMAN(sqrthslp_,nh_,yh_)
      eht_ = Gxr_*yh_*vh_

   end function CREHT

   !-----------------------------------------------------------------------

   function GRMR(wr_, rlen_) result(mr_)
      ! Function for estimate wetted channel area

      ! Arguments
      real(r8), intent(in) :: wr_   ! storage of water [m3]
      real(r8), intent(in) :: rlen_ ! channel length [m]
      real(r8)             :: mr_   ! wetted channel area [m2]

      mr_ = wr_ / rlen_
   end function GRMR

   !-----------------------------------------------------------------------

   function GRHT(mt_, twid_) result(ht_)
      ! Function for estimating water depth assuming rectangular channel

      ! Arguments
      real(r8), intent(in) :: mt_, twid_      ! wetted channel area [m2], channel width [m]
      real(r8)             :: ht_             ! water depth [m]

      if(mt_ <= TINYVALUE) then
         ht_ = 0._r8
      else
         ht_ = mt_ / twid_
      end if
   end function GRHT

   !-----------------------------------------------------------------------

   function GRPT(ht_, twid_) result(pt_)
      ! Function for estimating wetted perimeter assuming rectangular channel

      ! Arguments
      real(r8), intent(in) :: ht_, twid_      ! water depth [m], channel width [m]
      real(r8)             :: pt_             ! wetted perimeter [m]

      if(ht_ <= TINYVALUE) then
         pt_ = 0._r8
      else
         pt_ = twid_ + 2._r8 * ht_
      end if
   end function GRPT

   !-----------------------------------------------------------------------

   function GRRR(mr_, pr_) result(rr_)
      ! Function for estimating hydraulic radius

      ! Arguments
      real(r8), intent(in) :: mr_, pr_        ! wetted area [m2] and perimeter [m]
      real(r8)             :: rr_             ! hydraulic radius [m]

      if(pr_ <= TINYVALUE) then
         rr_ = 0._r8
      else
         rr_ = mr_ / pr_
      end if
   end function GRRR

   !-----------------------------------------------------------------------

   function GRHR(mr_, rwidth_, rwidth0_, rdepth_) result(hr_)
      ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
      ! here assuming the channel cross-section consists of three parts, from bottom to up,
      ! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
      ! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
      ! part 3 is a rectagular with the width rwid0

      ! Arguments
      real(r8), intent(in) :: mr_, rwidth_, rwidth0_, rdepth_ ! wetted channel area [m2], channel width [m], flood plain width [m], water depth [m]
      real(r8)             :: hr_                             ! water depth [m]

      ! Local variables
      real(r8) :: SLOPE1  ! slope of flood plain, TO DO
      real(r8) :: deltamr_

      SLOPE1 = SLOPE1def
      if(mr_ <= TINYVALUE) then
         hr_ = 0._r8
      else
         if(mr_ - rdepth_*rwidth_ <= TINYVALUE) then ! not flooded
            hr_ = mr_/rwidth_
         else ! if flooded, the find out the equivalent depth
            if(mr_ > rdepth_*rwidth_ + (rwidth_ + rwidth0_)*SLOPE1*((rwidth0_-rwidth_)/2._r8)/2._r8 + TINYVALUE) then
               deltamr_ = mr_ - rdepth_*rwidth_ - (rwidth_ + rwidth0_)*SLOPE1*((rwidth0_ - rwidth_)/2._r8)/2._r8;
               hr_ = rdepth_ + SLOPE1*((rwidth0_ - rwidth_)/2._r8) + deltamr_/(rwidth0_);
            else
               deltamr_ = mr_ - rdepth_*rwidth_;
               hr_ = rdepth_ + (-rwidth_+sqrt((rwidth_*rwidth_)+4._r8*deltamr_/SLOPE1))*SLOPE1/2._r8
            end if
         end if
      end if
   end function GRHR

   !-----------------------------------------------------------------------

   function GRPR(hr_, rwidth_, rwidth0_,rdepth_) result(pr_)
      ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
      ! here assuming the channel cross-section consists of three parts, from bottom to up,
      ! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
      ! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
      ! part 3 is a rectagular with the width rwid0

      ! Arguments
      real(r8), intent(in) :: hr_, rwidth_, rwidth0_, rdepth_ ! water depth [m], channel width [m], flood plain width [m], water depth [m]
      real(r8)             :: pr_                             ! water depth [m]

      ! Local variables
      real(r8) :: SLOPE1  ! slope of flood plain, TO DO
      real(r8) :: deltahr_
      logical, save :: first_call = .true.

      SLOPE1 = SLOPE1def
      if (first_call) then
         sinatanSLOPE1defr = 1.0_r8/(sin(atan(SLOPE1def)))
      endif
      first_call = .false.

      if(hr_ < TINYVALUE) then
         pr_ = 0._r8
      else
         if(hr_ <= rdepth_ + TINYVALUE) then ! not flooded
            pr_ = rwidth_ + 2._r8*hr_
         else
            if(hr_ > rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1 + TINYVALUE) then
               deltahr_ = hr_ - rdepth_ - ((rwidth0_-rwidth_)/2._r8)*SLOPE1
               pr_ = rwidth_ + 2._r8*(rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1*sinatanSLOPE1defr + deltahr_)
            else
               pr_ = rwidth_ + 2._r8*(rdepth_ + (hr_ - rdepth_)*sinatanSLOPE1defr)
            end if
         end if
      end if
   end function GRPR

   !-----------------------------------------------------------------------

   subroutine mosart_physics_restart()

      integer :: nt,nr

      associate( &
           ! ctl
           begr        => ctl%begr,            & ! local index start
           endr        => ctl%endr,            & ! local index end
           ntracers    => ctl%ntracers_tot,    & ! total number of tracers
           ! hillslope states
           wh          => TRunoff%wh,          & ! storage of surface water, [m]
           yh          => TRunoff%yh,          & ! depth of surface water, [m]
           ! subnetwork channel states
           wt          => Trunoff%wt,          & ! storage of surface water, [m3]
           yt          => Trunoff%yt,          & ! water depth, [m]
           mt          => Trunoff%mt,          & ! cross section area, [m2]
           pt          => Trunoff%pt,          & ! wetness perimeter, [m]
           rt          => TRunoff%rt,          & ! hydraulic radii, [m]
           ! main channel states
           wr          => Trunoff%wr,          & ! storage of surface water, [m3]
           mr          => Trunoff%mr,          & ! cross section area, [m2]
           yr          => Trunoff%yr,          & ! water depth. [m]
           pr          => Trunoff%pr,          & ! wetness perimeter, [m]
           rr          => Trunoff%rr           & ! hydraulic radius, [m]
      )

      do nt = 1,ntracers
         do nr = begr,endr
            call UpdateState_hillslope(wh(nr,nt), &          ! input
                 yh(nr,nt))                                  ! output
            call UpdateState_subnetwork(nr, wt(nr,nt), &     ! input
                 mt(nr,nt), yt(nr,nt), pt(nr,nt), rt(nr,nt)) ! output
            call UpdateState_mainchannel(nr, wr(nr,nt), &    ! input
                 mr(nr,nt), yr(nr,nt), pr(nr,nt), rr(nr,nt)) ! output
            ctl%volr(nr,nt) = (wt(nr,nt) + wr(nr,nt) + wh(nr,nt)*ctl%area(nr))
         enddo
      enddo

      end associate

   end subroutine mosart_physics_restart


end module mosart_physics
