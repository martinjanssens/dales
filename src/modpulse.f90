!> \file modpulse.f90
!!  Adds a pulse of qt to a single, horizontal, spatial wavenumber, 
!!  over a predefined height and at a preset time

!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

module modpulse

use modglobal, only : longint
implicit none
private
public :: initpulse, pulse, lpulse

save
! pulse governing variables (from namelist)
  logical :: lpulse     = .false.
  integer(kind=longint) :: timepulse
  real    :: amppulse_qt, radius_qt, amppulse_thl, radius_thl, amppulse_w, radius_w, zminpulse, zmaxpulse
  real, allocatable :: qtav0(:),qtav1(:), thlav0(:),thlav1(:), wmav1(:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initpulse
    use mpi
    use modmpi,    only :myid,D_MPI_BCAST, &
                         mpierr,commwrld
    use modglobal, only :cexpnr,runtime,ifnamopt,fname_options, &
                         checknamelisterror,tres,k1

    implicit none

    integer ierr
    namelist/NAMPULSE/ &
    lpulse, timepulse, amppulse_qt, radius_qt, amppulse_thl, radius_thl, amppulse_w, radius_w, zminpulse, zmaxpulse

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMPULSE,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMPULSE')
      write(6 ,NAMPULSE)
      close(ifnamopt)
    end if

    timepulse = timepulse/tres

    call D_MPI_BCAST(lpulse       ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(timepulse          ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(amppulse_qt        ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(radius_qt          ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(amppulse_thl       ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(radius_thl         ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(amppulse_w         ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(radius_w           ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(zminpulse          ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(zmaxpulse          ,1,0,commwrld,mpierr)
    
    allocate(qtav0(k1), qtav1(k1), thlav0(k1),thlav1(k1), wmav1(k1))

  end subroutine initpulse

  subroutine pulse

    use modglobal, only : rk3step,timee,dt_lim
    use modmpi,    only : myid
    implicit none
    if (.not. lpulse) return
    if (rk3step/=3) return

    if(timee<timepulse) then
      dt_lim = minval((/dt_lim,timepulse-timee/))
      return
    end if
    if (timee>=timepulse) then
      if (myid==0) then
        print *, 'Performing scalar pulse...'
      end if
      call do_pulse
      lpulse = .false.
    end if

  end subroutine pulse

  subroutine do_pulse

    use modfields     , only : qtm,thlm,wm
    use modglobal     , only : i1,j1,k1,imax,jmax, &
                               itot,jtot,dx,dy,zf,zh,pi, &
                               ih,jh,ijtot
    use modmpi        , only : myidx,myidy,myid,slabsum

    logical kstartflag
    integer kstart, kend, i, j, k
    real    xf, yf, qtpulse, thlpulse, wpulse, center_x, center_y

    kstartflag = .true.
    xf         = myidx*imax*dx
    yf         = myidy*jmax*dy
    qtpulse    = 0.
    thlpulse   = 0.
    wpulse     = 0.
    qtav0      = 0.
    qtav1      = 0.
    thlav0     = 0.
    thlav1     = 0.
    wmav1      = 0.

    ! Calculate the levels to apply the perturbation at
    do k=1,k1
      if (zf(k) >= zminpulse .and. kstartflag) then
        kstart = k
        kstartflag = .false.
      end if
      if (zf(k) > zmaxpulse) then
        kend = k
        exit
      end if
    end do
    if (myid == 0) then
      print *, 'kstart, kend', kstart, kend 
    end if

    ! Apply the perturbation - following protocol for cpmip

    ! Calculate domain-mean profiles    
    call slabsum(qtav0 ,1,k1,qtm ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(thlav0 ,1,k1,thlm ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qtav0 = qtav0 / ijtot
    thlav0 = thlav0 / ijtot

    ! Find the domain center
    center_x = dx*itot*0.5
    center_y = dy*jtot*0.5

    ! Add the pulse
    do i=2,i1
        xf = myidx*imax*dx + dx*(i - 1.5)
        do j=2,j1
        yf = myidy*jmax*dy + dy*(j - 1.5)
        
        ! qt pulse
        if ((xf-center_x)**2.0D0+(yf-center_y)**2.0D0 .lt. radius_qt**2.0D0) then
            qtpulse = amppulse_qt * cos( pi/2.0D0 * (sqrt( (xf-center_x)**2.0D0 + (yf-center_y)**2.0D0 ) / radius_qt) )**2.0D0
        else
            qtpulse = 0.
        end if

        ! thl pulse
        if ((xf-center_x)**2.0D0+(yf-center_y)**2.0D0 .lt. radius_thl**2.0D0) then
            thlpulse = amppulse_thl * cos( pi/2.0D0 * (sqrt( (xf-center_x)**2.0D0 + (yf-center_y)**2.0D0 ) / radius_thl) )**2.0D0
        else
            thlpulse = 0.
        end if

        ! w pulse
        if ((xf-center_x)**2.0D0+(yf-center_y)**2.0D0 .lt. radius_w**2.0D0) then
            wpulse = amppulse_w * cos( pi/2.0D0 * (sqrt( (xf-center_x)**2.0D0 + (yf-center_y)**2.0D0 ) / radius_w) )**2.0D0
        else
            wpulse = 0.
        end if

        ! if (myid == 0) then
        !     print *, 'x, y, qtpulse', xf, yf, qtpulse
        ! end if

        do k=kstart,kend 
            qtm(i,j,k) = qtm(i,j,k) + qtpulse
            thlm(i,j,k) = thlm(i,j,k) + thlpulse

            ! Apply in linearly increasing fashion, such that the (scaled) divergence
            ! amppulse_w/(z(kend)-z(kstart) is constant with height
            wm(i,j,k) = wm(i,j,k) + wpulse * (zh(k) - zh(kstart)) / (zh(kend) - zh(kstart))
        end do
        end do
    end do

    ! Calculate domain-mean profiles again
    call slabsum(qtav1  ,1,k1,qtm  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(thlav1 ,1,k1,thlm ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(wmav1  ,1,k1,wm   ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qtav1 = qtav1 / ijtot
    thlav1 = thlav1 / ijtot
    wmav1 = wmav1 / ijtot

    ! And subtract the difference everywhere
    do i=2,i1
        do j=2,j1
        do k=kstart,kend 
            qtm(i,j,k) = qtm(i,j,k) - (qtav1(k) - qtav0(k))
            thlm(i,j,k) = thlm(i,j,k) - (thlav1(k) - thlav0(k))
            wm(i,j,k) = wm(i,j,k) - wmav1(k)
        end do
        end do
    end do

  end subroutine do_pulse

end module modpulse