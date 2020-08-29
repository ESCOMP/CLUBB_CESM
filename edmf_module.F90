module edmf_module

  use clubb_precision, only: core_rknd

  implicit none
  private
  save

  public :: integrate_mf

  contains

  ! =============================================================================== !
  !  Eddy-diffusivity mass-flux routine                                             !
  ! =============================================================================== !

  subroutine integrate_mf( nz, zt, dzt, p, iexner, nup, u, v,                       & ! input
                           th, thl, thl_zm, thv, qt, qt_zm,                         & ! input
                           ust, wthl, wqt, pblh, qc,                                & ! input
                           dry_a,   moist_a,                                        & ! output - plume diagnostics
                           dry_w,   moist_w,                                        & ! output - plume diagnostics
                           dry_qt,  moist_qt,                                       & ! output - plume diagnostics
                           dry_thl, moist_thl,                                      & ! output - plume diagnostics
                           dry_u,   moist_u,                                        & ! output - plume diagnostics
                           dry_v,   moist_v,                                        & ! output - plume diagnostics
                                    moist_qc,                                       & ! output - variables needed for solver
                           ae, aw, awthl, awqt, awql, awqi, awu, awv,               & ! output - diagnosed fluxes BEFORE mean field update
                           thlflx, qtflx,                                           &
                           edmf_enti )

  ! Original author: Marcin Kurowski, JPL
  ! Modified heavily by Mikael Witte, UCLA/JPL for implementation in CESM2/E3SM
  ! Additional Modifications by Adam Herrington, NCAR

  !
  ! Variables needed for solver
  ! ae = sum_i (1-a_i)
  ! aw3 = sum (a_i w_i)
  ! aws3 = sum (a_i w_i*s_i); s=thl*cp
  ! aws3,awqv3,awql3,awqi3,awu3,awv3 similar as above except for different variables
  !

  ! - mass flux variables are computed on edges (i.e. momentum grid):
  !  upa,upw,upqt,... kts:kte+1
  !  dry_a,moist_a,dry_w,moist_w, ... kts:kte+1

     use physconst,          only: rair, cpair, gravit, latvap, latice, zvir
     use parameters_tunable, only: &
         mf_L0, &
         mf_ent0, & 
         mf_wa, &
         mf_wb   

     integer, intent(in) :: nz,nup
     real(kind=core_rknd), dimension(nz), intent(in) :: u, v, th, thl, qt, qc, thv, zt, dzt, iexner, p ! thermodynamic/midpoint levels
     real(kind=core_rknd), dimension(nz), intent(in) :: thl_zm, qt_zm
     real(kind=core_rknd), intent(in)                :: ust,wthl,wqt
     real(kind=core_rknd), intent(inout)             :: pblh

     real(kind=core_rknd),dimension(nz), intent(out) :: dry_a, moist_a, dry_w, moist_w, &
                                                        dry_qt, moist_qt, dry_thl, moist_thl, &
                                                        dry_u, moist_u, dry_v, moist_v, moist_qc
     real(kind=core_rknd),dimension(nz), intent(out) :: ae, aw, awthl, awqt, awql, awqi, awu, awv
     real(kind=core_rknd),dimension(nz), intent(out) :: thlflx, qtflx
     real(kind=core_rknd),dimension(nz), intent(out) :: edmf_enti     

     ! INTERNAL VARIABLES

     ! sums over all plumes
     real(kind=core_rknd), dimension(nz)     :: moist_th, dry_th, awqv, awth

     ! updraft properties
     real(kind=core_rknd), dimension(nz,nup) :: upw, upthl, upqt, upqc, upth, upqv, upql,  &
                                                upqi, upa, upu, upv, upthv, ups
     ! entrainment variables
     real(kind=core_rknd), dimension(nz,nup) :: ent,entf
     integer,  dimension(nz,nup)             :: enti

     ! other variables
     integer                                 :: k,i,ih
     real(kind=core_rknd)                    :: wthv, wstar, qstar, thstar, sigmaw, sigmaqt, sigmath, z0, &
                                                wmin, wmax, wlv, wtv, wp
     real(kind=core_rknd)                    :: b, qtn, thln, thvn, thn, qcn, qln, qin, un, vn, wn2, &
                                                entexp, entexpu, entw

     ! internal surface cont
     real(kind=core_rknd)                    :: iexh

     ! parameters defining initial conditions for updrafts
     real(kind=core_rknd),parameter          :: &
                                                pwmin = 1.5,&
                                                pwmax = 3.

     ! min values to avoid singularities
     real(kind=core_rknd),parameter          :: &
                                                wstarmin = 1.e-3, &
                                                pblhmin  = 100.

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!! BEGIN CODE !!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! INITIALIZE OUTPUT VARIABLES
     ! set updraft properties to zero
     dry_a     = 0.
     moist_a   = 0.
     dry_w     = 0.
     moist_w   = 0.
     dry_qt    = 0.
     moist_qt  = 0.
     dry_thl   = 0.
     moist_thl = 0.
     dry_u     = 0.
     moist_u   = 0.
     dry_v     = 0.
     moist_v   = 0.
     moist_qc  = 0.
     ! outputs - variables needed for solver
     aw        = 0.
     ! aws       = 0.
     awthl     = 0.
     awqt      = 0.
     awqv      = 0.
     awql      = 0.
     awqi      = 0.
     awu       = 0.
     awv       = 0.
     thlflx    = 0.
     qtflx     = 0.

     ! this is the environmental area - by default 1.
     ae = 1.

     ! START MAIN COMPUTATION
     upw   = 0.
     upthl = 0.
     upthv = 0.
     upqt  = 0.
     upa   = 0.
     upu   = 0.
     upv   = 0.
     upqc  = 0.
     ent   = 0.
     upth  = 0.
     upql  = 0.
     upqi  = 0.
     upqv  = 0.

     pblh = max(pblh,pblhmin)
     wthv = wthl+zvir*thv(1)*wqt

     ! if surface buoyancy is positive then do mass-flux, otherwise not
     if ( wthv > 0.0 ) then

       ! compute entrainment coefficient
       ! get dz/L0
       do i=1,nup
         do k=2,nz
           entf(k,i) = dzt(k) / mf_L0
         enddo
       enddo

       ! get Poisson P(dz/L0)
       call poisson( nz, nup, entf, enti, thl(nz))

       edmf_enti(:) = enti(:,1)

       ! entrainent: Ent=Ent0/dz*P(dz/L0)
       do i=1,nup
         do k=2,nz
           ent(k,i) = real( enti(k,i))*mf_ent0/dzt(k) ! MKW TODO: also input invrs_dzt? only used here.
         enddo
       enddo

       ! surface conditions
       wstar  = max( wstarmin, (gravit/thv(1)*wthv*pblh)**(1./3.) ) ! MKW NOTE: is it better to use kts instead of index 1 for surface?
       qstar  = wqt / wstar
       thstar = wthl / wstar

       sigmaw  = 0.572 * wstar     / 1.
       sigmaqt = 2.89 * abs(qstar) / 1.
       sigmath = 2.89 * abs(thstar)/ 1.

       wmin = sigmaw * pwmin
       wmax = sigmaw * pwmax

       do i=1,nup

         wlv = wmin + (wmax-wmin) / (real(nup)) * (real(i)-1.)
         wtv = wmin + (wmax-wmin) / (real(nup)) * real(i)

         upw(1,i) = 0.5 * (wlv+wtv)
         upa(1,i) = 0.5 * erf( wtv/(sqrt(2.)*sigmaw) ) &
                    - 0.5 * erf( wlv/(sqrt(2.)*sigmaw) )

         upu(1, i) = u(1)
         upv(1, i) = v(1)

         upqc(1,i)  = 0.
         upqt(1,i)  = qt(1)  + 0.32 * upw(1,i) * sigmaqt/sigmaw
         upthv(1,i) = thv(1) + 0.58 * upw(1,i) * sigmath/sigmaw
         upthl(1,i) = upthv(1,i) / (1.+zvir*upqt(1,i))
         upth(1,i)  = upthl(1,i)
         upqv(1,i)  = upqt(1,i)

       enddo

       ! integrate updrafts
       do i=1,nup
         do k=2,nz!kts+1,kte+1

           entexp  = exp(-ent(k,i)*dzt(k))
           entexpu = exp(-ent(k,i)*dzt(k)/3.)

           qtn  = qt(k-1) *(1.-entexp ) + upqt (k-1,i)*entexp
           thln = thl(k-1)*(1.-entexp ) + upthl(k-1,i)*entexp
           un   = u(k-1)  *(1.-entexpu) + upu  (k-1,i)*entexpu
           vn   = v(k-1)  *(1.-entexpu) + upv  (k-1,i)*entexpu

           iexh = (1.e5 / p(k))**(rair/cpair) ! MKW NOTE: why not just use CLUBB exner??
           call condensation_mf(qtn, thln, p(k), iexh, &
                                 thvn, qcn, thn, qln, qin)

           b=gravit*(0.5*(thvn+upthv(k-1,i))/thv(k-1)-1.)
           !b=mapl_grav*(thvn/thv(k-1)-1.)

           ! Wn^2
           ! to avoid singularities w equation has to be computed diferently if wp==0
           wp = mf_wb*ent(k,i)
           if (wp==0.) then
             wn2 = upw(k-1,i)**2+2.*mf_wa*b*dzt(k-1)
           else
             entw = exp(-2.*wp*dzt(k-1))
             wn2 = entw*upw(k-1,i)**2+mf_wa*b/(mf_wb*ent(k,i))*(1.-entw)
           end if

           if (wn2>0.) then
             upw(k,i)   = sqrt(wn2)
             upthv(k,i) = thvn
             upthl(k,i) = thln
             upqt(k,i)  = qtn
             upqc(k,i)  = qcn
             upu(k,i)   = un
             upv(k,i)   = vn
             upa(k,i)   = upa(k-1,i)
             upth(k,i)  = thn
             upql(k,i)  = qln
             upqi(k,i)  = qin
             upqv(k,i)  = qtn - qcn
           else
             exit
           end if
         enddo
       enddo

       ! writing updraft properties for output
       ! all variables, except Areas are now multipled by the area
       do k=1,nz

         ! first sum over all i-updrafts
         do i=1,nup
           if (upqc(k,i)>0.) then
             moist_a(k)   = moist_a(k)   + upa(k,i)
             moist_w(k)   = moist_w(k)   + upa(k,i)*upw(k,i)
             moist_qt(k)  = moist_qt(k)  + upa(k,i)*upqt(k,i)
             moist_thl(k) = moist_thl(k) + upa(k,i)*upthl(k,i)
             moist_th(k)  = moist_th(k)  + upa(k,i)*upth(k,i)
             moist_u(k)   = moist_u(k)   + upa(k,i)*upu(k,i)
             moist_v(k)   = moist_v(k)   + upa(k,i)*upv(k,i)
             moist_qc(k)  = moist_qc(k)  + upa(k,i)*upqc(k,i)
           else
             dry_a(k)     = dry_a(k)     + upa(k,i)
             dry_w(k)     = dry_w(k)     + upa(k,i)*upw(k,i)
             dry_qt(k)    = dry_qt(k)    + upa(k,i)*upqt(k,i)
             dry_thl(k)   = dry_thl(k)   + upa(k,i)*upthl(k,i)
             dry_th(k)    = dry_th(k)    + upa(k,i)*upth(k,i)
             dry_u(k)     = dry_u(k)     + upa(k,i)*upu(k,i)
             dry_v(k)     = dry_v(k)     + upa(k,i)*upv(k,i)
           endif
         enddo

         if ( dry_a(k) > 0. ) then
           dry_w(k)   = dry_w(k)   / dry_a(k)
           dry_qt(k)  = dry_qt(k)  / dry_a(k)
           dry_thl(k) = dry_thl(k) / dry_a(k)
           dry_th(k)  = dry_th(k)  / dry_a(k)
           dry_u(k)   = dry_u(k)   / dry_a(k)
           dry_v(k)   = dry_v(k)   / dry_a(k)
         else
           dry_w(k)   = 0.
           dry_qt(k)  = 0.
           dry_thl(k) = 0.
           dry_th(k)  = 0.
           dry_u(k)   = 0.
           dry_v(k)   = 0.
         endif

         if ( moist_a(k) > 0. ) then
           moist_w(k)   = moist_w(k)   / moist_a(k)
           moist_qt(k)  = moist_qt(k)  / moist_a(k)
           moist_thl(k) = moist_thl(k) / moist_a(k)
           moist_th(k)  = moist_th(k)  / moist_a(k)
           moist_u(k)   = moist_u(k)   / moist_a(k)
           moist_v(k)   = moist_v(k)   / moist_a(k)
           moist_qc(k)  = moist_qc(k)  / moist_a(k)
         else
           moist_w(k)   = 0.
           moist_qt(k)  = 0.
           moist_thl(k) = 0.
           moist_th(k)  = 0.
           moist_u(k)   = 0.
           moist_v(k)   = 0.
           moist_qc(k)  = 0.
         endif

       enddo

       do k=1,nz
         do i=1,nup
           ae  (k) = ae  (k) - upa(k,i)
           aw  (k) = aw  (k) + upa(k,i)*upw(k,i)
           awu (k) = awu (k) + upa(k,i)*upw(k,i)*upu(k,i)
           awv (k) = awv (k) + upa(k,i)*upw(k,i)*upv(k,i)
           !aws (k) = aws (k) + upa(k,i)*upw(k,i)*upth(k,i)*cpair
           !aws (k) = aws (k) + upa(k,i)*upw(k,i)*ups(k,i)
           awthl(k)= awthl(k)+ upa(k,i)*upw(k,i)*upthl(k,i) !*cpair/iexh
           awth(k) = awth(k) + upa(k,i)*upw(k,i)*upth(k,i) !*cpair/iexh
           awqt(k) = awqt(k) + upa(k,i)*upw(k,i)*upqt(k,i)
           awqv(k) = awqv(k) + upa(k,i)*upw(k,i)*upqv(k,i)
           awql(k) = awql(k) + upa(k,i)*upw(k,i)*upql(k,i)
           awqi(k) = awqi(k) + upa(k,i)*upw(k,i)*upqi(k,i)
         enddo
       enddo

       do k=2,nz
         ! iexh = (1.e5/p(k))**(rair/cpair)
         thlflx(k)= awthl(k) - aw(k)*thl_zm(k) ! MKW NOTE: used to be slflx, but CLUBB works on thl
         !sflx( k)= (awth(k) - aw(k)*0.5*(th(k-1)+th(k)) )*cpair/iexh ! not using this since all s/sl stuff is handled in clubb_cam_tend
         qtflx(k)= awqt(k)  - aw(k)*qt_zm(k)
       enddo
       !iexh = (1.e5/p(kts))**(rair/cpair)
       thlflx(1) = 0.
       !sflx(kts)  = 0.
       qtflx(1) = 0.

       print*,'max(1-ae)=',maxval(1.-ae)

     end if  ! ( wthv > 0.0 )

  end subroutine integrate_mf

  subroutine condensation_mf( qt, thl, p, iex, thv, qc, th, ql, qi)
  !
  ! zero or one condensation for edmf: calculates thv and qc
  !
     use physconst,       only: cpair, zvir
     use wv_saturation,      only : qsat

     real(kind=core_rknd),intent(in) :: qt,thl,p,iex
     real(kind=core_rknd),intent(out):: thv,qc,th,ql,qi

     !local variables
     integer :: niter,i
     real(kind=core_rknd) :: diff,t,qs,qcold,es,wf

     ! max number of iterations
     niter=50
     ! minimum difference
     diff=2.e-5

     qc=0.
     t=thl/iex

     !by definition:
     ! T   = Th*Exner, Exner=(p/p0)^(R/cp)   (1)
     ! Thl = Th - L/cp*ql/Exner              (2)
     !so:
     ! Th  = Thl + L/cp*ql/Exner             (3)
     ! T   = Th*Exner=(Thl+L/cp*ql/Exner)*Exner    (4)
     !     = Thl*Exner + L/cp*ql
     do i=1,niter
       wf = get_watf(t)
       t = thl/iex+get_alhl(wf)/cpair*qc   !as in (4)

       ! qsat, p is in pascal (check!)
       call qsat(t,p,es,qs)
       qcold = qc
       qc = max(0.5*qc+0.5*(qt-qs),0.)
       if (abs(qc-qcold)<diff) exit
     enddo

     wf = get_watf(t)
     t = thl/iex+get_alhl(wf)/cpair*qc


     call qsat(t,p,es,qs)
     qc = max(qt-qs,0.)
     thv = (thl+get_alhl(wf)/cpair*iex*qc)*(1.+zvir*(qt-qc)-qc)
     th = t*iex
     qi = qc*(1.-wf)
     ql = qc*wf

     contains

     function get_watf(t)
       real(kind=core_rknd) :: t,get_watf,tc
       real(kind=core_rknd), parameter :: &
       tmax=-10., &
       tmin=-40.

       tc=t-273.16

       if (tc>tmax) then
         get_watf=1.
       else if (tc<tmin) then
         get_watf=0.
       else
         get_watf=(tc-tmin)/(tmax-tmin);
       end if

     end function get_watf


     function get_alhl(wf)
     !latent heat of the mixture based on water fraction
       use physconst,        only : latvap , latice
       real(kind=core_rknd) :: get_alhl,wf

       get_alhl = wf*latvap+(1.-wf)*(latvap+latice)

     end function get_alhl

  end subroutine condensation_mf

  subroutine poisson(nz,nup,lambda,poi,state)

       integer, intent(in) :: nz,nup
       real(kind=core_rknd), intent(in) :: state
       real(kind=core_rknd), dimension(nz,nup), intent(in) :: lambda
       integer, dimension(nz,nup), intent(out) :: poi
       integer :: i,j

       call set_seed_from_state(state)

       do i=1,nz
         do j=1,nup
           call knuth(lambda(i,j),poi(i,j))
         enddo
       enddo

  end subroutine poisson

  subroutine set_seed_from_state(state)

       implicit none
       real(kind=core_rknd),intent(in) :: state
       integer, allocatable :: seed(:)
       integer :: i,n,tmpseed

       call random_seed(size = n)

       if (allocated(seed)) deallocate(seed)
       allocate(seed(n))

       tmpseed = int((state - int(state)) * 1000000000)
       do i=1,n
         seed(i) = tmpseed
       end do

       call random_seed(put=seed)
       deallocate(seed)

  end subroutine set_seed_from_state

  subroutine knuth(lambda,kout)
  !**********************************************************************
  ! Discrete random poisson from Knuth 
  ! The Art of Computer Programming, v3, 137-138
  !**********************************************************************
       implicit none

       real(kind=core_rknd), intent(in) :: lambda
       integer, intent(out) :: kout

       !Local variables
       real(kind=core_rknd) :: puni, tmpuni, explam
       integer :: k

       k = 0
       explam = exp(-1.*lambda)
       puni = 1.0
       do while (puni.gt.explam)
         k = k + 1
         call random_number(tmpuni)
         puni = puni*tmpuni
       end do
       kout = k - 1

  end subroutine knuth

end module edmf_module
