#include "rundeck_opts.h"

      subroutine obio_update(vrbos,kmax,i,j)
 
c  Performs updating of biological particles, uses mid-point
c  leap frog method.

      USE obio_dim
      USE obio_com, only: P_tend,obio_deltat,D_tend,C_tend
     .                   ,obio_P,det,car,dp1d,p1d
#ifdef TRACERS_Alkalinity
     .                   ,A_tend,alk1d
#ifdef TOPAZ_params
     .                   ,Ca_tend,ca_det_calc1d
#endif
#endif
#ifdef TRACERS_Ocean_O2
     .                   ,O_tend,o21d
#endif

      implicit none

      integer :: i,j,k
 
      integer :: nt,kmax
      real    :: Pnew,Dnew,Cnew,Anew,Canew,O2new
      logical :: vrbos
 
c  Loop to update
c   indexes are mixed because new H has already been computed
c   in update.F, but P has not been updated yet

      do 1000 k = 1,kmax

        do nt = 1,ntyp
         Pnew = (obio_P(k ,nt) +  P_tend(k,nt)*obio_deltat)
         obio_P(k,nt) = max(0.d0,pnew)
        enddo
 
        do nt = 1,ndet
         Dnew = (det(k,nt)     +  D_tend(k,nt)*obio_deltat)
          det(k,nt) = max(0.d0,Dnew)
        enddo

        do nt = 1,ncar
         Cnew = (car(k,nt) +  C_tend(k,nt)*obio_deltat)
         car(k,nt) = max(0.d0,Cnew) 
        enddo

#ifdef TRACERS_Alkalinity
         nt=1
         Anew = (alk1d(k) +  A_tend(k)*obio_deltat)
         alk1d(k) = max(0.d0,Anew)
#ifdef TOPAZ_params
         nt = 2
         Canew = (ca_det_calc1d(k) +  Ca_tend(k)*obio_deltat)
         ca_det_calc1d(k) = Canew
#endif
#endif

#ifdef TRACERS_Ocean_O2
         O2new = (o21d(k) +  O_tend(k)*obio_deltat)
         o21d(k) = O2new
#endif

 1000 continue


      return
      end
