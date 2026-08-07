! ###########################################################################################
!> \file atmos_coupling.F90
!> Procedures for coupling the MPAS dynamical core to the CCPP Physics.
!>
! ###########################################################################################
module atmos_coupling_mod
  use mpas_kind_types, only : mpas_kind => RKIND
  use ufs_mpas_io,     only : domain_ptr
  
  implicit none
  public :: ufs_physics_to_mpas
  public :: ufs_mpas_to_physics
  public :: ufs_microphysics_to_mpas
  public :: ufs_mpas_to_microphysics
  public :: ufs_mpas_grid_to_physics
  public :: ufs_mpas_sfc_to_physics
  
contains
  !> #########################################################################################
  !> Procedure to populate CCPP data containers with MPAS pool data.
  !> Called BEFORE CCPP Radiation and Physics Groups.
  !>
  !> Analogous to MPAS_to_physics in src/core_atmosphere/physics/mpas_atmphys_interface.F
  !>
  !> This procedure accesses MPAS data using MPAS native procedures and stores the data
  !> locally in the data-containers defined above. The MPAS "state" is then translated to the
  !> CCPP "state" needed by the physics.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_to_physics(physics_state, surface_state)
    use GFS_typedefs,         only : GFS_statein_type, GFS_sfcprop_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use atm_core,             only : atm_compute_output_diagnostics
    use mpas_kind_types,      only : RKIND
    use mpas_constants,       only : gravity, rvord

    ! Arguments
    type(GFS_statein_type),   intent(inout) :: physics_state
    type(GFS_sfcprop_type),   intent(inout) :: surface_state

    ! Locals
    type(mpas_pool_type), pointer :: state_pool, diag_pool, mesh_pool
    integer :: iCol, iLay, iTracer, ithread
    integer, pointer :: nCellsSolve, num_scalars, nVertLevels, index_qv
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(kind=RKIND) :: rho1, rho2, tem1, tem2, theta, fzm_p, fzp_p, z0, z1, z2, w1, w2,rho_a
    real(kind=RKIND), pointer :: qv(:,:), qc(:,:), qr(:,:), qi(:,:), qs(:,:), qg(:,:)
    real(kind=RKIND), pointer :: ux(:,:), uy(:,:), theta_m(:,:), rho_zz(:,:), zgrid(:,:), zz(:,:)
    real(kind=RKIND), pointer :: exner(:,:), tracers(:,:,:), pressure_b(:,:), pressure_p(:,:)
    real(kind=RKIND), pointer :: w(:,:), surface_pressure(:),  prsi(:,:), prsl(:,:), rho(:,:), dz(:,:)
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_to_physics'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)

    ! Get MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(diag_pool,  'uReconstructZonal',      ux)
    call mpas_pool_get_array(diag_pool,  'uReconstructMeridional', uy)
    call mpas_pool_get_array(state_pool, 'scalars',                tracers, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'w',                      w, timeLevel=1)
    call mpas_pool_get_array(diag_pool,  'exner',                  exner)
    call mpas_pool_get_array(mesh_pool,  'zgrid',                  zgrid)
    call mpas_pool_get_array(mesh_pool,  'zz',                     zz)
    call mpas_pool_get_array(state_pool, 'theta_m',                theta_m, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'rho_zz',                 rho_zz,  timeLevel=1)
    call mpas_pool_get_array(diag_pool,  'pressure_base',          pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p',             pressure_p)
    call mpas_pool_get_array(diag_pool,  'surface_pressure',       surface_pressure)

    ! Local variables
    allocate(prsl(nCellsSolve, nVertLevels))
    allocate(prsi(nCellsSolve, nVertLevels + 1))
    allocate(rho( nCellsSolve, nVertLevels))
    allocate(dz(  nCellsSolve, nVertLevels))
    
    ! Copy fields from MPAS data containers to physics data containers.
    ! [k, i] -> [i, k]
    ! Retain bottom-up convention
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels
             ! Scalars (tracer,layer,col) -> (col,layer,tracer)
             do iTracer = 1,num_scalars
                physics_state % qgrs(iCol,iLay,iTracer) = max(0._RKIND, tracers(iTracer,iLay,iCol))
             end do

             ! Air denisty (rho) (TODO: Pass to CCPP Physics)
             rho(iCol,iLay) = zz(iLay,iCol) * rho_zz(iLay,iCol)

             ! Potential temperature (theta_m -> theta)
             theta = theta_m(iLay,iCol) / (1._RKIND + rvord * tracers(index_qv,iLay,iCol))

             ! Air temperature (theta -> temp)
             physics_state % tgrs(iCol,iLay)   = theta*exner(iLay,iCol)

             ! Winds at grid center
             physics_state % ugrs(iCol,iLay)   = ux(iLay,iCol)
             physics_state % vgrs(iCol,iLay)   = uy(iLay,iCol)

             ! Layer-height
             physics_state % phil(iCol,iLay)   = 0.5*(zgrid(iLay+1,iCol)+zgrid(iLay,iCol))

             ! Level height
             physics_state % phii(iCol,iLay)   = zgrid(iLay,iCol)

             ! Layer thickness (TODO: Pass to CCPP physics)
             dz(iCol,iLay) = zgrid(iLay+1,iCol) - zgrid(iLay,iCol)

             ! Exner funciton
             physics_state % prslk(iCol,iLay)  = exner(iLay,iCol)

             ! MPAS provides vertical velocity at interfaces, compute layer mean.
             physics_state % vvl(iCol,iLay) = 0.5*(w(iLay,iCol) + w(iLay+1,iCol))

             ! Pressure (non-hydrostatic)
             prsl(iCol,iLay) = pressure_p(iLay,iCol) + pressure_b(iLay,iCol)

          end do
          ! Set surface temperature to lowest level temperature (revisit for coupling)
          theta = theta_m(1,iCol) / (1._RKIND + rvord * tracers(index_qv,1,iCol))
          surface_state % tsfc(iCol) = theta*exner(1,iCol)
       end do
    end do

    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          tem1 = zgrid(2,iCol) - zgrid(1,iCol)
          tem2 = zgrid(3,iCol) - zgrid(2,iCol)
          rho1 = rho_zz(1,iCol) * zz(1,iCol) * (1. + tracers(index_qv,1,iCol))
          rho2 = rho_zz(2,iCol) * zz(2,iCol) * (1. + tracers(index_qv,2,iCol))
          surface_pressure(iCol) = 0.5*gravity*(zgrid(2,iCol) - zgrid(1,iCol)) &
               * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
          surface_pressure(iCol) = surface_pressure(iCol) + pressure_p(1,iCol) + pressure_b(1,iCol)
       end do
    end do

    ! Interpolation of pressure and temperature from layer-center to layer-interface
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 2,nVertLevels
             tem1 = 1./(zgrid(iLay+1,iCol)-zgrid(iLay-1,iCol))
             fzm_p = (zgrid(iLay,  iCol)-zgrid(iLay-1,iCol)) * tem1
             fzp_p = (zgrid(iLay+1,iCol)-zgrid(iLay,  iCol)) * tem1
             physics_state % tgri(iCol,iLay) = fzm_p*physics_state % tgrs(iCol,iLay) + fzp_p*physics_state % tgrs(iCol,iLay-1)
             prsi(iCol,iLay) = fzm_p*prsl(iCol,iLay) + fzp_p*prsl(iCol,iLay-1)
          enddo
       enddo
    enddo

    ! Interpolation of pressure and temperature to the top-of-the-model
    iLay = nVertLevels + 1
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          z0 = zgrid(iLay,iCol)
          z1 = 0.5*(zgrid(iLay  ,iCol)+zgrid(iLay-1,iCol))
          z2 = 0.5*(zgrid(iLay-1,iCol)+zgrid(iLay-2,iCol))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.-w1
          physics_state % tgri(iCol,iLay) = w1*physics_state % tgrs(iCol,iLay-1) + w2*physics_state % tgrs(iCol,iLay-2)
          prsi(iCol,iLay) = exp(w1*log(prsl(iCol,iLay-1))+w2*log(prsl(iCol,iLay-2)))
       end do
    end do

    ! Recalculate the pressure and temperature  at the surface as an extrapolation of
    ! the pressures in the 2 layers above the surface.
    iLay = 1
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          z0 = zgrid(iLay,iCol)
          z1 = 0.5*(zgrid(iLay  ,iCol)+zgrid(iLay+1,iCol))
          z2 = 0.5*(zgrid(iLay+1,iCol)+zgrid(iLay+2,iCol))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.-w1
          physics_state % tgri(iCol,iLay) = w1*physics_state % tgrs(iCol,iLay) + w2*physics_state % tgrs(iCol,iLay+1)
          prsi(iCol,iLay) = w1*prsl(iCol,iLay)+w2*prsl(iCol,iLay+1)
       end do
    end do

    ! Calculation of the hydrostatic pressure
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          ! Pressure at layer-interfaces
          iLay = nVertLevels + 1
          physics_state % prsi(iCol,iLay) = prsi(iCol,iLay)
          do iLay = nVertLevels,1,-1
             rho_a = rho(iCol,iLay) / (1.+tracers(index_qv,iLay,iCol))
             physics_state % prsi(iCol,iLay)  = physics_state % prsi(iCol,iLay+1) + &
                  gravity*rho(iCol,iLay)*dz(iCol,iLay)
          end do
          ! Pressure at layer-centers
          do iLay = nVertLevels,1,-1
             physics_state % prsl(iCol,iLay) = 0.5*(physics_state % prsi(iCol,iLay+1) + physics_state % prsi(iCol,iLay) )
          end do
          ! surface pressure:
          physics_state % pgr(iCol) = physics_state % prsi(iCol,1)
       end do
    end do

    ! Housekeeping
    deallocate (prsl)
    deallocate (prsi)
    deallocate (dz)
    deallocate (rho)
    nullify (mesh_pool)
    nullify (state_pool)
    nullify (diag_pool)

  end subroutine ufs_mpas_to_physics

  !> #########################################################################################
  !> Procedure to update MPAS state with physics (CCPP) tendencies.
  !> Called AFTER physics, BEFORE calling dynamics (current timestep).
  !>
  !> - Convert from theta to theta_m (tend_theta_phys)
  !> - Update dynamic tendencies. (tend_theta_dyn)
  !> - Update scalar tendencies (tend_scalars_dyn)
  !>
  !> Analogous to phys_get_tend in physics/mpas_atmphys_todynamics.F
  !> Here, instead of updating the state with physics tendencies from the MPAS "tend_pool", we
  !> will use tendencies from the CCPP Physics data containers.
  !>
  !> #########################################################################################
  subroutine ufs_physics_to_mpas(physics_state)
    use GFS_typedefs,       only : GFS_stateout_type
    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use mpas_kind_types,    only : RKIND
    use mpas_constants,     only : rv, rgas, gravity

    ! Arguments
    type(GFS_stateout_type), intent(in) :: physics_state

    ! Locals
    type(mpas_pool_type),  pointer :: state_pool, mesh_pool, tend_pool, diag_pool
    real(kind=RKIND), pointer :: mass(:,:), exner(:,:), theta_m(:,:), zgrid(:,:), zz(:,:)
    real(kind=RKIND), pointer :: pressure_b(:,:), pressure_p(:,:), tend_th_phys(:,:)
    real(kind=RKIND), pointer :: tend_theta_phys(:,:), tend_theta_dyn(:,:)
    real(kind=RKIND), pointer :: scalars(:,:,:), tend_scalars_phys(:,:,:), tend_scalars_dyn(:,:,:)
    real(kind=RKIND), pointer :: surface_pressure(:)
    integer, pointer :: nCellsSolve, num_scalars, nVertLevels 
    integer, pointer :: index_qv => null()
    integer, pointer :: index_qc => null()
    integer, pointer :: index_qi => null()
    integer, pointer :: index_qr => null()
    integer, pointer :: index_qs => null()
    integer, pointer :: index_qg => null()
    integer, pointer :: index_nc => null()
    integer, pointer :: index_ni => null()
    integer, pointer :: index_nifa => null()
    integer, pointer :: index_nwfa => null()
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    integer :: iCol,iLay,ithread,iScalar
    real(kind=RKIND):: coeff, tem1, tem2, rho1, rho2
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_physics_to_mpas'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend',  tend_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)

    ! Get MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)
    call mpas_pool_get_dimension(state_pool, 'index_qc',    index_qc)
    call mpas_pool_get_dimension(state_pool, 'index_qi',    index_qi)
    call mpas_pool_get_dimension(state_pool, 'index_qr',    index_qr)
    call mpas_pool_get_dimension(state_pool, 'index_qs',    index_qs)
    call mpas_pool_get_dimension(state_pool, 'index_qg',    index_qg)
    call mpas_pool_get_dimension(state_pool, 'index_nc',    index_nc)
    call mpas_pool_get_dimension(state_pool, 'index_ni',    index_ni)
    call mpas_pool_get_dimension(state_pool, 'index_nifa',  index_nifa)
    call mpas_pool_get_dimension(state_pool, 'index_nwfa',  index_nwfa)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(state_pool,'theta_m',          theta_m,1)
    call mpas_pool_get_array(state_pool,'scalars',          scalars,1)
    call mpas_pool_get_array(state_pool,'rho_zz',           mass,   1)
    call mpas_pool_get_array(mesh_pool, 'zz',               zz)
    call mpas_pool_get_array(mesh_pool, 'zgrid',            zgrid)
    call mpas_pool_get_array(diag_pool, 'surface_pressure', surface_pressure)
    call mpas_pool_get_array(diag_pool, 'pressure_base',    pressure_b)
    call mpas_pool_get_array(diag_pool, 'pressure_p',       pressure_p)
    call mpas_pool_get_array(diag_pool, 'exner',            exner)

    ! Allocate local variables
    allocate(tend_th_phys(nVertLevels,nCellsSolve+1))
    allocate(tend_theta_phys(nVertLevels,nCellsSolve+1))
    allocate(tend_scalars_phys(num_scalars, nVertLevels,nCellsSolve+1))
    tend_th_phys(:,:)        = 0._RKIND
    tend_theta_phys(:,:)     = 0._RKIND
    tend_scalars_phys(:,:,:) = 0._RKIND

    ! GJF: Add accumulated tendencies from the physics group
    do ithread=1,nThreads
      do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
        do iLay = 1,nVertLevels 
          tend_th_phys(iLay,iCol) = tend_th_phys(iLay,iCol) + (physics_state % dtdt(iCol,iLay)/exner(iLay,iCol))*mass(iLay,iCol)
          tend_scalars_phys(index_qv,iLay,iCol) = tend_scalars_phys(index_qv,iLay,iCol) + physics_state % dqdt(iCol,iLay,index_qv)*mass(iLay,iCol)
        end do
      end do
    end do

    if(associated(index_qc)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
            tend_scalars_phys(index_qc,iLay,iCol) = tend_scalars_phys(index_qc,iLay,iCol) + physics_state % dqdt(iCol,iLay,index_qc)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    if(associated(index_qi)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
            tend_scalars_phys(index_qi,iLay,iCol) = tend_scalars_phys(index_qi,iLay,iCol) + physics_state % dqdt(iCol,iLay,index_qi)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    if(associated(index_qr)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
            tend_scalars_phys(index_qr,iLay,iCol) = tend_scalars_phys(index_qr,iLay,iCol) + physics_state % dqdt(iCol,iLay,index_qr)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    if(associated(index_qs)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
            tend_scalars_phys(index_qs,iLay,iCol) = tend_scalars_phys(index_qs,iLay,iCol) + physics_state % dqdt(iCol,iLay,index_qs)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    if(associated(index_qg)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
            tend_scalars_phys(index_qg,iLay,iCol) = tend_scalars_phys(index_qg,iLay,iCol) + physics_state % dqdt(iCol,iLay,index_qg)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    if(associated(index_nc)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
            tend_scalars_phys(index_nc,iLay,iCol) = tend_scalars_phys(index_nc,iLay,iCol) + physics_state % dqdt(iCol,iLay,index_nc)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    if(associated(index_ni)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
            tend_scalars_phys(index_ni,iLay,iCol) = tend_scalars_phys(index_ni,iLay,iCol) + physics_state % dqdt(iCol,iLay,index_ni)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    if(associated(index_nifa)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
            tend_scalars_phys(index_nifa,iLay,iCol) = tend_scalars_phys(index_nifa,iLay,iCol) + physics_state % dqdt(iCol,iLay,index_nifa)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    if(associated(index_nwfa)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
            tend_scalars_phys(index_nwfa,iLay,iCol) = tend_scalars_phys(index_nwfa,iLay,iCol) + physics_state % dqdt(iCol,iLay,index_nwfa)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    ! Convert from potential temperature to modified potential temperature (theta -> theta_m)
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1, nVertLevels
             coeff = (1. + rv/rgas * scalars(index_qv,iLay,iCol))
             tend_th_phys(iLay,iCol) = coeff * tend_th_phys(iLay,iCol) + rv/rgas * theta_m(iLay,iCol) * tend_scalars_phys(index_qv,iLay,iCol) / coeff
             tend_theta_phys(iLay,iCol) = tend_theta_phys(iLay,iCol) + tend_th_phys(iLay,iCol)
          end do
       end do
    end do

    ! Update MPAS state tendencies
    call mpas_pool_get_array(tend_pool, 'theta_m', tend_theta_dyn)
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1, nVertLevels
               tend_theta_dyn(iLay,iCol) = tend_theta_dyn(iLay,iCol) + tend_theta_phys(iLay,iCol)
          end do
       end do
    end do

    call mpas_pool_get_array(tend_pool, 'scalars_tend', tend_scalars_dyn)
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1, nVertLevels
             do iScalar = 1, num_scalars
               tend_scalars_dyn(iScalar,iLay,iCol) = tend_scalars_dyn(iScalar,iLay,iCol) + tend_scalars_phys(iScalar,iLay,iCol)
             end do
          end do
       end do
    end do

    ! Housekeeping
    deallocate(tend_th_phys)
    deallocate(tend_theta_phys)
    deallocate(tend_scalars_phys)
    
    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          tem1 = zgrid(2,iCol) - zgrid(1,iCol)
          tem2 = zgrid(3,iCol) - zgrid(2,iCol)
          rho1 = mass(1,iCol) * zz(1,iCol) * (1. + scalars(index_qv,1,iCol))
          rho2 = mass(2,iCol) * zz(2,iCol) * (1. + scalars(index_qv,2,iCol))
          surface_pressure(iCol) = 0.5*gravity*( zgrid(2,iCol) -  zgrid(1,iCol)) &
               * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
          surface_pressure(iCol) = surface_pressure(iCol) + pressure_p(1,iCol) + &
               pressure_b(1,iCol)
       enddo
    enddo

    ! Housekeeping
    nullify (state_pool)
    nullify (mesh_pool)
    nullify (diag_pool)
  end subroutine ufs_physics_to_mpas

  !> #########################################################################################
  !> Procedure to compute diabatic heating tendency from microphysics and store in MPAS pool.
  !> Called AFTER microphysics, BEFORE calling dynamics (next timestep)
  !>
  !> This routine is where we compute the diabatic heating rate, "rt_diabatic_tend", due to the
  !> microphysics.
  !> Follows microphysics_to_MPAS from src/core_atmosphere/physics/mpas_atmphys_interface.F.
  !> compute the diabatic heating rate, "rt_diabatic_tend" and save in MPAS memory for use by
  !> the dynamics at the subsequent time step in dynamics/mpas_atm_time_integration.F.
  !> Additionally, update any other fields needed by the dynamics (e.g., theta_m, rtheta_p)
  !>
  !> #########################################################################################
  subroutine ufs_microphysics_to_mpas(physics_state)
    use GFS_typedefs,       only : GFS_stateout_type
    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array
    use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_config
    use mpas_constants,     only : gravity, rvord, rv, rgas, p0, cp
    use mpas_kind_types,    only : RKIND

    ! Arguments
    type(GFS_stateout_type),     intent(in   ) :: physics_state

    ! Locals
    type(mpas_pool_type), pointer :: diag_pool, mesh_pool, state_pool, tend_pool
    integer :: iCol, ithread, iLay, iTracer
    integer, pointer :: nCellsSolve
    integer, pointer :: index_qv => null()
    integer, pointer :: index_qc => null()
    integer, pointer :: index_qi => null()
    integer, pointer :: index_qr => null()
    integer, pointer :: index_qs => null()
    integer, pointer :: index_qg => null()
    integer, pointer :: index_nc => null()
    integer, pointer :: index_ni => null()
    integer, pointer :: index_nr => null()
    integer, pointer :: index_ns => null()
    integer, pointer :: index_ng => null()
    integer, pointer :: index_nifa => null()
    integer, pointer :: index_nwfa => null()
    integer, pointer :: num_scalars, nVertLevels
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(kind=RKIND) :: rho1, rho2, tem1, tem2, coeff, rcv, theta_dyn
    real(kind=RKIND), pointer :: config_dt
    real(kind=RKIND), pointer :: tracers(:,:,:), rt_diabatic_tend(:,:), rho_zz(:,:), theta_m(:,:)
    real(kind=RKIND), pointer :: zz(:,:), zgrid(:,:), exner(:,:), exner_b(:,:), rtheta_b(:,:), theta(:,:)
    real(kind=RKIND), pointer :: rtheta_p(:,:), pressure_b(:,:), pressure_p(:,:), surface_pressure(:), dtheta_dt_mp(:,:)
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_microphysics_to_mpas'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend',  tend_pool)

    ! Model timestep
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_dt', config_dt)

    ! MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)
    call mpas_pool_get_dimension(state_pool, 'index_qc',    index_qc)
    call mpas_pool_get_dimension(state_pool, 'index_qi',    index_qi)
    call mpas_pool_get_dimension(state_pool, 'index_qr',    index_qr)
    call mpas_pool_get_dimension(state_pool, 'index_qs',    index_qs)
    call mpas_pool_get_dimension(state_pool, 'index_qg',    index_qg)
    call mpas_pool_get_dimension(state_pool, 'index_nc',    index_nc)
    call mpas_pool_get_dimension(state_pool, 'index_ni',    index_ni)
    call mpas_pool_get_dimension(state_pool, 'index_nr',    index_nr)
    call mpas_pool_get_dimension(state_pool, 'index_ns',    index_ns)
    call mpas_pool_get_dimension(state_pool, 'index_ng',    index_ng)
    call mpas_pool_get_dimension(state_pool, 'index_nifa',  index_nifa)
    call mpas_pool_get_dimension(state_pool, 'index_nwfa',  index_nwfa)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(state_pool, 'scalars',          tracers, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'rho_zz',           rho_zz,  timeLevel=1)
    call mpas_pool_get_array(state_pool, 'theta_m',          theta_m, timeLevel=1)
    call mpas_pool_get_array(mesh_pool,  'zgrid',            zgrid)
    call mpas_pool_get_array(mesh_pool,  'zz',               zz)
    call mpas_pool_get_array(diag_pool,  'pressure_base',    pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p',       pressure_p)
    call mpas_pool_get_array(diag_pool,  'surface_pressure', surface_pressure)
    call mpas_pool_get_array(diag_pool,  'exner',            exner)
    call mpas_pool_get_array(diag_pool,  'exner_base',       exner_b)
    call mpas_pool_get_array(diag_pool,  'rtheta_p',         rtheta_p)
    call mpas_pool_get_array(diag_pool,  'rtheta_base',      rtheta_b)
    call mpas_pool_get_array(tend_pool,  'rt_diabatic_tend', rt_diabatic_tend)
    call mpas_pool_get_array(diag_pool,  'dtheta_dt_mp',     dtheta_dt_mp)

    ! The MPAS version of microphysics schemes update the state within the dynamics;
    ! for CCPP/UFS, we will need to update the state variables here for use by the dynamics.
    ! Also, Update water vapor, cloud liquid water, rain mixing ratios, modified potential temperature,
    ! and potential temperature heating rate from microphysics
    allocate(theta(nVertLevels, nCellsSolve))
    rcv = rgas/(cp-rgas)
    do ithread=1,nThreads
      do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
        do iLay = 1,nVertLevels
          
          ! Initialize diabatic heating tendency (initially fill with theta_m before updating)
          rt_diabatic_tend(iLay,iCol) = theta_m(iLay,iCol)
          
          ! Update potential temperature (theta) with microphysics tendency
          coeff = (1._RKIND + rvord * tracers(index_qv,iLay,iCol))
          theta_dyn = theta_m(ilay,iCol)/coeff
          theta(iLay,iCol) = theta_dyn + config_dt * (physics_state % dtdt(iCol,iLay) / exner(iLay,iCol))
          
          ! Scalars (col,layer,tracer) -> (tracer,layer,col)
          do iTracer = 1,num_scalars
            tracers(iTracer,iLay,iCol) = max(0._RKIND, tracers(iTracer,iLay,iCol) + config_dt * physics_state % dqdt(iCol,iLay,iTracer))
          end do
          
          ! update the virtual temperature coefficient with updated qv
          coeff = (1._RKIND + rvord * tracers(index_qv,iLay,iCol))
          ! Modified potential temperature (theta ->theta_m)
          theta_m(iLay,iCol) = theta(iLay,iCol)*coeff
          
          ! Now compute diabatic heating due to microphsyics, save for next time step
          rt_diabatic_tend(iLay,iCol) = (theta_m(iLay,iCol) - rt_diabatic_tend(iLay,iCol)) / config_dt
          
          ! Save the straight theta tendency due to microphysics
          dtheta_dt_mp(iLay,iCol) =  (theta(iLay,iCol) - theta_dyn) / config_dt
          
          ! Density weighted perturbation potential temperature
          rtheta_p(iLay,iCol) = rho_zz(iLay,iCol) * theta_m(iLay,iCol) - rtheta_b(iLay,iCol)
          
          ! Exner function
          exner(iLay,iCol) = (zz(iLay,iCol)*(rgas/P0)*(rtheta_p(iLay,iCol)+rtheta_b(iLay,iCol)))**rcv

          ! Perturbation pressure
          pressure_p(iLay,iCol) = zz(iLay,iCol)*rgas*(exner(iLay,iCol)*rtheta_p(iLay,iCol) + &
                                    (exner(iLay,iCol)-exner_b(iLay,iCol))*rtheta_b(iLay,iCol))  
        end do
      end do
    end do
    
    ! write(*,*) 'num_scalars',num_scalars
    ! if (associated(index_qv)) write(*,*) 'mean max/min ten qv',sum(tracers(index_qv,:,:)) / real(size(tracers(index_qv,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qv)), config_dt*minval(physics_state % ten_q(:,:,index_qv))
    ! if (associated(index_qc)) write(*,*) 'mean max/min ten qc',sum(tracers(index_qc,:,:)) / real(size(tracers(index_qc,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qc)), config_dt*minval(physics_state % ten_q(:,:,index_qc)) 
    ! if (associated(index_qi)) write(*,*) 'mean max/min ten qi',sum(tracers(index_qi,:,:)) / real(size(tracers(index_qi,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qi)), config_dt*minval(physics_state % ten_q(:,:,index_qi))
    ! if (associated(index_qr)) write(*,*) 'mean max/min ten qr',sum(tracers(index_qr,:,:)) / real(size(tracers(index_qr,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qr)), config_dt* minval(physics_state % ten_q(:,:,index_qr))
    ! if (associated(index_qs)) write(*,*) 'mean max/min ten qs',sum(tracers(index_qs,:,:)) / real(size(tracers(index_qs,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qs)), config_dt*minval(physics_state % ten_q(:,:,index_qs))
    ! if (associated(index_qg)) write(*,*) 'mean max/min ten qg',sum(tracers(index_qg,:,:)) / real(size(tracers(index_qg,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qg)), config_dt*minval(physics_state % ten_q(:,:,index_qg))
    ! if (associated(index_nc)) write(*,*) 'mean max/min ten nc',sum(tracers(index_nc,:,:)) / real(size(tracers(index_nc,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_nc)), config_dt*minval(physics_state % ten_q(:,:,index_nc))
    ! if (associated(index_ni)) write(*,*) 'mean max/min ten ni',sum(tracers(index_ni,:,:)) / real(size(tracers(index_ni,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_ni)), config_dt*minval(physics_state % ten_q(:,:,index_ni))
    ! if (associated(index_nr)) write(*,*) 'mean max/min ten nr',sum(tracers(index_nr,:,:)) / real(size(tracers(index_nr,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_nr)), config_dt*minval(physics_state % ten_q(:,:,index_nr))
    ! if (associated(index_ns)) write(*,*) 'mean max/min ten ns',sum(tracers(index_ns,:,:)) / real(size(tracers(index_ns,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_ns)), config_dt*minval(physics_state % ten_q(:,:,index_ns))
    ! if (associated(index_ng)) write(*,*) 'mean max/min ten ng',sum(tracers(index_ng,:,:)) / real(size(tracers(index_ng,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_ng)), config_dt*minval(physics_state % ten_q(:,:,index_ng))
    ! if (associated(index_nifa)) write(*,*) 'mean max/min ten nifa',sum(tracers(index_nifa,:,:)) / real(size(tracers(index_nifa,:,:))), config_dt*maxval(physics_state % ten_q(:,:,index_nifa)), config_dt*minval(physics_state % ten_q(:,:,index_nifa))
    ! if (associated(index_nwfa)) write(*,*) 'mean max/min ten nwfa',sum(tracers(index_nwfa,:,:)) /real(size(tracers(index_nwfa,:,:))),config_dt*maxval(physics_state % ten_q(:,:,index_nwfa)), config_dt*minval(physics_state % ten_q(:,:,index_nwfa))
    

    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          tem1 = zgrid(2,iCol) - zgrid(1,iCol)
          tem2 = zgrid(3,iCol) - zgrid(2,iCol)
          rho1 = rho_zz(1,iCol) * zz(1,iCol) * (1. + tracers(index_qv,1,iCol))
          rho2 = rho_zz(2,iCol) * zz(2,iCol) * (1. + tracers(index_qv,2,iCol))
          surface_pressure(iCol) = 0.5*gravity*(zgrid(2,iCol) - zgrid(1,iCol)) &
               * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
          surface_pressure(iCol) = surface_pressure(iCol) + pressure_p(1,iCol) + pressure_b(1,iCol)
       end do
    end do

    ! Housekeeping
    nullify (state_pool)
    nullify (mesh_pool)
    nullify (diag_pool)

  end subroutine ufs_microphysics_to_mpas

  !> #########################################################################################
  !> Procedure to update physics (CCPP) state using updated MPAS state.
  !> Called AFTER dynamics, BEFORE microphysics.
  !> 
  !> Analogous to microphysics_from_MPAS in src/core_atmosphere/physics/mpas_atmphys_interface.F
  !>
  !> #########################################################################################
  subroutine ufs_mpas_to_microphysics(physics_state, physics_statein)
    use GFS_typedefs,         only : GFS_stateout_type, GFS_statein_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use mpas_constants,       only : rvord
    use mpas_kind_types,    only : RKIND

    ! Arguments
    type(GFS_stateout_type), intent(inout) :: physics_state
    type(GFS_statein_type),  intent(inout) :: physics_statein

    ! Locals
    type(mpas_pool_type), pointer :: state_pool, diag_pool, mesh_pool
    integer :: ithread, iCol, iTracer, iLay
    integer, pointer :: num_scalars, nVertLevels, nCellsSolve, index_qv, index_qc, index_qr, index_qi, index_qs, index_qg, index_ni, index_nr
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(kind=RKIND), pointer :: rho_zz(:,:), theta_m(:,:), zz(:,:), zgrid(:,:), exner(:,:)
    real(kind=RKIND), pointer :: tracers(:,:,:), rho(:,:), pressure_b(:,:), pressure_p(:,:), w(:,:)
    real(kind=RKIND) :: theta, pres, z, dz
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_to_microphysics'

    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)

    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)
    call mpas_pool_get_dimension(state_pool, 'index_qc',    index_qc)
    call mpas_pool_get_dimension(state_pool, 'index_qr',    index_qr)
    call mpas_pool_get_dimension(state_pool, 'index_qi',    index_qi)
    call mpas_pool_get_dimension(state_pool, 'index_qs',    index_qs)
    call mpas_pool_get_dimension(state_pool, 'index_qg',    index_qg)
    call mpas_pool_get_dimension(state_pool, 'index_ni',    index_ni)
    call mpas_pool_get_dimension(state_pool, 'index_nr',    index_nr)
    

    call mpas_pool_get_array(state_pool, 'rho_zz',       rho_zz,  timeLevel=1)
    call mpas_pool_get_array(state_pool, 'theta_m',      theta_m, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'scalars',      tracers, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'w',            w,       timeLevel=1)
    call mpas_pool_get_array(mesh_pool,  'zz',           zz)
    call mpas_pool_get_array(mesh_pool,  'zgrid',        zgrid)
    call mpas_pool_get_array(diag_pool,  'exner',        exner)
    call mpas_pool_get_array(diag_pool,  'pressure_base',pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p'   ,pressure_p)

    allocate(rho(nCellsSolve, nVertLevels))
    ! Update fields needed by microphysics...
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels
             ! Scalars (tracer,layer,col) -> (col,layer,tracer)
             do iTracer = 1,num_scalars
                physics_state % gq0(iCol,iLay,iTracer) = max(0._RKIND, tracers(iTracer,iLay,iCol))
             end do

             ! Air denisty (rho) (TODO: Pass to CCPP Physics)
             rho(iCol,iLay) = zz(iLay,iCol) * rho_zz(iLay,iCol)

             ! Potential temperature (theta_m -> theta)
             theta = theta_m(iLay,iCol) / (1._RKIND + rvord * max(0._RKIND,tracers(index_qv,iLay,iCol)))

             ! Air temperature (theta -> t)
             physics_state % gt0(iCol,iLay) = theta*exner(iLay,iCol)

             ! Pressure
             pres = pressure_b(iLay,iCol) + pressure_p(iLay,iCol)

             ! Height and layer-thickness (TODO: Pass to CCPP Physics)
             z  = zgrid(iLay,iCol)
             dz = zgrid(iLay+1,iCol) - zgrid(iLay,iCol)

             ! Vertical velocity
             physics_statein % vvl(iCol,iLay) = 0.5*(w(iLay,iCol) + w(iLay+1,iCol))
          end do
       end do
    end do

    deallocate(rho)
    nullify(diag_pool)
    nullify(mesh_pool)
    nullify(state_pool)

    ! GJF: Remove microphysics heating from state before calling microphysics. This is done
    ! at line 3317 of mpas_atm_time_integration.F/atm_recover_large_step_variables_work.
 
  end subroutine ufs_mpas_to_microphysics

!> #########################################################################################
!> Procedure to transfer MPAS grid information to physics DDTs.
!>
!> #########################################################################################
  subroutine ufs_mpas_grid_to_physics(physics_grid)
    use GFS_typedefs,         only : GFS_grid_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension
    use mpas_pool_routines,   only : mpas_pool_get_array, mpas_pool_get_config
    use mpas_kind_types,      only : RKIND
    use mpas_constants,       only : pii
    use mpas_log,             only : mpas_log_write
    use mpas_derived_types,   only : MPAS_LOG_ERR, MPAS_LOG_WARN, MPAS_LOG_CRIT
    ! Arguments
    type(GFS_grid_type),      intent(inout) :: physics_grid

    ! Locals
    type(mpas_pool_type), pointer :: mesh_pool
    integer :: i, ierr, ithread
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(RKIND), pointer :: lat(:), lon(:), area(:), meshDensity(:)
    real(RKIND), pointer :: nominalMinDc, config_len_disp
    real(RKIND)          :: rad2deg
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_grid_to_physics'

    ierr = 0
    rad2deg = 180.0_RKIND/pii
    
    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    
    call mpas_pool_get_array(mesh_pool,  'latCell',     lat)
    call mpas_pool_get_array(mesh_pool,  'lonCell',     lon)
    call mpas_pool_get_array(mesh_pool,  'areaCell',    area)
    call mpas_pool_get_array(mesh_pool,  'meshDensity', meshDensity)
    
    ! (from mpas_atm_core.F/atm_core_init Determine horizontal length scale used by horizontal diffusion and 3-d divergence damping
    nullify(nominalMinDc)
    call mpas_pool_get_array(mesh_pool, 'nominalMinDc', nominalMinDc)

    nullify(config_len_disp)
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_len_disp', config_len_disp)

    ! If config_len_disp was specified as a valid value, use that
    if (config_len_disp > 0.0_RKIND) then
      ! But if nominalMinDc was available in the input file and is different, print a warning
      if (nominalMinDc > 0.0_RKIND .and. abs(nominalMinDc - config_len_disp) > 1.0e-6_RKIND * config_len_disp) then
        call mpas_log_write(subname // ' WARNING: nominalMinDc was read from input file as a positive value ($r) that differs', &
                            realArgs=[nominalMinDc], messageType=MPAS_LOG_WARN)
        call mpas_log_write(subname // ' WARNING: from the specified config_len_disp value ($r)', &
                            realArgs=[config_len_disp], messageType=MPAS_LOG_WARN)
      end if
      nominalMinDc = config_len_disp
    ! Otherwise, try to use nominalMinDc
    else
      if (nominalMinDc > 0.0_RKIND) then
        call mpas_log_write('Setting config_len_disp to $r based on nominalMinDc value in input file', realArgs=[nominalMinDc])
          config_len_disp = nominalMinDc
      else
         call mpas_log_write(subname // ' ERROR: Both config_len_disp and nominalMinDc are <= 0.0.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write(subname // ' ERROR: Please either specify config_len_disp in the &nhyd_model namelist group,', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write(subname // ' ERROR: or use an input file that provides a valid value for the nominalMinDc variable.', &
                             messageType=MPAS_LOG_ERR)
        ierr = 1
      end if
    end if
    if (ierr/=0)  call mpas_log_write(subname // ' ERROR: Call to ufs_mpas_grid_to_physics() failed', messageType=MPAS_LOG_CRIT)

    do ithread = 1,nThreads
       do i = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          physics_grid % xlat(i)   = lat(i)
          physics_grid % xlon(i)   = lon(i)
          physics_grid % xlat_d(i) = physics_grid % xlat(i) * rad2deg
          physics_grid % xlon_d(i) = physics_grid % xlon(i) * rad2deg
         ! write(*,*) 'ithread, i, xlat, xlon, xlat_d, xlon_d',ithread, i, physics_grid % xlat(i),physics_grid % xlon(i),physics_grid % xlat_d(i),physics_grid % xlon_d(i)
          physics_grid % sinlat(i) = sin(physics_grid % xlat(i))
          physics_grid % coslat(i) = sqrt(1.0_RKIND - physics_grid % sinlat(i) * physics_grid % sinlat(i))
          physics_grid % area(i)   = area(i)
          !formula for dx comes from mpas_atmphys_driver_gwdo.F instead of sqrt(area) as in FV3
          physics_grid % dx(i)     = config_len_disp / meshDensity(i)**0.25
       end do
    end do
    
  end subroutine ufs_mpas_grid_to_physics

!> #########################################################################################
!> Procedure to transfer MPAS information to physics srfprop DDT
!>
!> #########################################################################################
  subroutine ufs_mpas_sfc_to_physics(physics_sfcprop)
    use GFS_typedefs,         only : GFS_sfcprop_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
    use mpas_kind_types,      only : RKIND

    ! Arguments
    type(GFS_sfcprop_type),      intent(inout) :: physics_sfcprop
    ! Locals
    type(mpas_pool_type), pointer :: mesh_pool, sfc_input
    integer :: i, ierr, iCol, ithread
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(RKIND), pointer :: landmask(:), sst(:), snow(:), tmn(:), albbck(:)
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_sfc_to_physics'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input', sfc_input)

    !using fv3atm_sfc_io.F90/Sfc_io_transfer() as a template; mpas_init_atm_static.F from MPAS-model for syntax
    call mpas_pool_get_array(mesh_pool, 'landmask',  landmask)
    call mpas_pool_get_array(sfc_input, 'sst',       sst)
    call mpas_pool_get_array(sfc_input, 'snow',      snow)
    call mpas_pool_get_array(sfc_input, 'tmn' ,      tmn)
    call mpas_pool_get_array(sfc_input, 'sfc_albbck',albbck)
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          !physics_sfcprop % slmsk(iCol) = landmask(iCol)
          physics_sfcprop % tsfco(iCol) = sst(iCol)
          physics_sfcprop % weasd(iCol) = snow(iCol)
          physics_sfcprop % tg3(iCol)   = tmn(iCol)
          !zorl - z0/znt in MPAS, read in as sfz0; landuse_init_forMPAS not used yet; diag_physics pool, z0 - not initialized yet
          !alvsf, alvwf, alnsf, alnwf - MPAS doesn't split into visible/nir and strong/weak coszen dependency; set these to the value that we have (background snow-free albedo of surface)?
          physics_sfcprop % alvsf(iCol) = albbck(iCol)
          physics_sfcprop % alvwf(iCol) = albbck(iCol)
          physics_sfcprop % alnsf(iCol) = albbck(iCol)
          physics_sfcprop % alnwf(iCol) = albbck(iCol)
       end do
    end do

  end subroutine ufs_mpas_sfc_to_physics

  !> #########################################################################################
  !> Procedure to populate MPAS diag_phys pool with CCPP data.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_phys_diag(radiation,diagnostics)
    use GFS_typedefs,         only : GFS_radtend_type
    use GFS_typedefs,         only : GFS_diag_type
    use mpas_kind_types,      only : RKIND
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_config

    ! Arguments
    type(GFS_radtend_type), intent(in) :: radiation
    type(GFS_diag_type),    intent(in) :: diagnostics

    ! Locals
    type(mpas_pool_type), pointer :: diag_phys
    real(RKIND),dimension(:),pointer :: swdnb,swdnbc,swupb,swupbc
    real(RKIND),dimension(:),pointer :: lwdnb,lwdnbc,lwupb,lwupbc
    real(RKIND),dimension(:,:),pointer :: refl10cm
    real(RKIND),dimension(:),pointer :: rainc,rainnc,frainnc,snownc,graupelnc
    real(RKIND),dimension(:),pointer :: raincv,rainncv,snowncv,graupelncv
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    integer :: iCol, ithread
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_phys_diag'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics',  diag_phys)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(diag_phys,'swdnb'     , swdnb     )
    call mpas_pool_get_array(diag_phys,'swdnbc'    , swdnbc    )
    call mpas_pool_get_array(diag_phys,'swupb'     , swupb     )
    call mpas_pool_get_array(diag_phys,'swupbc'    , swupbc    )
    call mpas_pool_get_array(diag_phys,'lwdnb'     , lwdnb     )
    call mpas_pool_get_array(diag_phys,'lwdnbc'    , lwdnbc    )
    call mpas_pool_get_array(diag_phys,'lwupb'     , lwupb     )
    call mpas_pool_get_array(diag_phys,'lwupbc'    , lwupbc    )
    call mpas_pool_get_array(diag_phys,'refl10cm'  , refl10cm  )
    call mpas_pool_get_array(diag_phys,'rainc'     , rainc     )
    call mpas_pool_get_array(diag_phys,'rainnc'    , rainnc    )
    call mpas_pool_get_array(diag_phys,'frainnc'   , frainnc   )
    call mpas_pool_get_array(diag_phys,'snownc'    , snownc    )
    call mpas_pool_get_array(diag_phys,'graupelnc' , graupelnc )
    call mpas_pool_get_array(diag_phys,'raincv'    , raincv    )
    call mpas_pool_get_array(diag_phys,'rainncv'   , rainncv   )
    call mpas_pool_get_array(diag_phys,'snowncv'   , snowncv   )
    call mpas_pool_get_array(diag_phys,'graupelncv', graupelncv)

    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          ! Radiation fluxes at surface
          swdnb(iCol)  = radiation%sfcfsw(iCol)%dnfxc
          swdnbc(iCol) = radiation%sfcfsw(iCol)%dnfx0
          swupb(iCol)  = radiation%sfcfsw(iCol)%upfxc
          swupbc(iCol) = radiation%sfcfsw(iCol)%upfx0
          lwdnb(iCol)  = radiation%sfcflw(iCol)%dnfxc
          lwdnbc(iCol) = radiation%sfcflw(iCol)%dnfx0
          lwupb(iCol)  = radiation%sfcflw(iCol)%upfxc
          lwupbc(iCol) = radiation%sfcflw(iCol)%upfx0
          ! Reflectivity
          refl10cm(:,iCol) = diagnostics%refl_10cm(iCol,:)
          ! Instantaneous precipitation
          raincv(iCol)     = diagnostics%rain(iCol)
          rainncv(iCol)    = diagnostics%rainc(iCol)
          snowncv(iCol)    = diagnostics%snow(iCol)
          graupelncv(iCol) = diagnostics%graupel(iCol)
          ! Accumulated precipitation
          rainc(iCol)      = diagnostics%cnvprcp(iCol)
          rainnc(iCol)     = diagnostics%totprcp(iCol)
          frainnc(iCol)    = diagnostics%totice(iCol)
          snownc(iCol)     = diagnostics%totsnw(iCol)
          graupelnc(iCol)  = diagnostics%totgrp(iCol)
       end do
    end do
  end subroutine ufs_mpas_phys_diag  
end module atmos_coupling_mod
