***
      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,vk,snstar,
     &                r2,fallback,sigmahold,kick_info,disrupt,bkick)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'

      integer kw,k,snstar,sn,safety
      real*8 m1,m2,m1n,mbi,mbf,mdif
      real*8 ecc,sep,sepn,jorb,ecc2
      real*8 pi,twopi,gmrkm,yearsc,rsunkm
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mm,em,dif,der,del,r
      real*8 u1,u2,vk,v(4),s,theta,phi
      real*8 sphi,cphi,stheta,ctheta,salpha,calpha
      real*8 x_tilt,y_tilt,z_tilt
      real*8 mu,cmu,smu,omega,comega,somega
      real*8 cmu1,smu1,comega1,somega1
      real*8 vr,vr2,vk2,vn2,hn2
      real*8 vs(3),v1,v2,v3
      real*8 mx1,mx2,r2
      real*8 sigmah,RotInvX
      real*8 signs,sigc,psins,psic,cpsins,spsins,cpsic,spsic
      real*8 csigns
      real*8 semilatrec,cangleofdeath,angleofdeath,energy
      real*8 fallback,sigmahold,bound
      real*8 mean_mns,mean_mej,alphakick,betakick
      real*8 bkick(20)
      real*8 maxwellian
      INTEGER directcollapse,ECS
      real*8 ffb,mfin
      COMMON /KICKSN/ ffb,directcollapse,ECS,mfin
* Output
      logical output,disrupt
*
      real*8 kick_info(2,17)
      real ran3,xx
      external ran3
*
      output = .false. !useful for debugging...
      safety = 0

      do k = 1,3
         vs(k) = 0.d0
      enddo
*     if(kw.eq.14.and.bhflag.eq.0) goto 95
*
* To take into account the ECS
      if(ECS .eq. 1)then
* Check if the supernova is driven by electron-cupture... 
          maxwellian = 15   
      else
* ... or iron core collapse. 
          maxwellian = 15
      endif
*
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
      gmrkm = 1.906125d+05
*
* Find the initial separation by randomly choosing a mean anomaly.
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
         xx = RAN3(idum1)
         mm = xx*twopi
         em = mm
 2       dif = em - ecc*SIN(em) - mm
         if(ABS(dif/mm).le.1.0d-04) goto 3
         der = 1.d0 - ecc*COS(em)
         del = dif/der
         em = em - del
         goto 2
 3       continue
         r = sep*(1.d0 - ecc*COS(em))
*
* Find the initial relative velocity vector.
         salpha = SQRT((sep*sep*(1.d0-ecc*ecc))/(r*(2.d0*sep-r)))
         calpha = (-1.d0*ecc*SIN(em))/SQRT(1.d0-ecc*ecc*COS(em)*COS(em))
         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
         vr = SQRT(vr2)
      else
         vr = 0.d0
         vr2 = 0.d0
         salpha = 0.d0
         calpha = 0.d0
      endif
*
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
* Use Henon''s method for pairwise components (Douglas Heggie 22/5/97).
      do 20 k = 1,2
         u1 = RAN3(idum1)
         u2 = RAN3(idum1)
* Generate two velocities from polar coordinates S & THETA.
* way to avoid negative values in LOG
         ss = 1.d0 - u1
         if(ss.le.0.d0)then
            ss = 1.d-6
         endif
         s = maxwellian*SQRT(-2.d0*LOG(ss))
         theta = twopi*u2
         v(2*k-1) = s*COS(theta)
         v(2*k) = s*SIN(theta)
 20   continue
      vk2 = v(1)**2 + v(2)**2 + v(3)**2
      vk = SQRT(vk2)
*
***
* Module kick velocity for Black Holes...
* Old case:
      if(kw.eq.14.and.bhflag.eq.1)then
         mrem = m1n
         vk = vk*(1.d0 - ffb)
         vk2 = vk*vk
*     ... if they are generated via direct collapse vk = 0
***
* module vick velocity for Black Holes and Neutron Stars
* Fallback case:
      elseif((kw.eq.14.or.kw.eq.13).and.(bhflag.eq.2))then
         mrem = m1n
         vk = vk*(1.d0 - ffb)
         vk2 = vk*vk
***
* Generate Kick Velocity for Black Holes and Neutron Stars
* Giacobbo&Mapelli 2020 case:
      elseif((kw.eq.14.or.kw.eq.13).and.(bhflag.eq.3))then
         mrem = m1n
*
         if(kw.eq.14 .and. ffb.eq.1.d0)then
            mfin = mrem
         endif  
* mean mass ejected <mej> in Msun (in a population of 10^6 single stars)
         if(remnantflag.eq.3) mean_ej = 9.d0
         if(remnantflag.eq.4) mean_ej = 7.5d0
         vk = vk*((mfin-mrem)/mean_ej)*(1.2d0/mrem)
         vk2 = vk*vk
***   
* NO-modulation of kick kelocity for Black Holes and Neutron Stars
* Full kick case:
      elseif((kw.eq.14.or.kw.eq.13).and.(bhflag.eq.4))then
         mrem = m1n
         vk = vk
         vk2 = vk*vk  
***
      elseif((kw.eq.14.and.bhflag.eq.0).or.kw.lt.0)then
         vk2 = 0.d0
         vk = 0.d0
         if(kw.lt.0) kw = 13
***     
      endif
***
* Correction (Noticed by Peter & Mirek)      
      u1 = RAN3(idum1)
      u2 = RAN3(idum1)
      theta = twopi*u2
***
      sphi = -1.d0 + 2.d0*u1
      phi = ASIN(sphi)
      cphi = COS(phi)
      stheta = SIN(theta)
      ctheta = COS(theta)
      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90
*
* Determine the magnitude of the new relative velocity.
* (NG 12/17: there was a wrong sign in the brackets)
      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha+stheta*cphi*calpha)
* Calculate the new semi-major axis.
      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
      sep = 1.d0/sep
*     if(sep.le.0.d0)then
*        ecc = 1.1d0
*        goto 90
*     endif
* Determine the magnitude of the cross product of the separation vector
* and the new relative velocity.
      v1 = vk2*sphi*sphi
      v2 = (vk*ctheta*cphi-vr*salpha)**2
      hn2 = r*r*(v1 + v2)
* Calculate the new eccentricity.
      ecc2 = 1.d0 - hn2/(gmrkm*sep*(m1n+m2))
      ecc2 = MAX(ecc2,0.d0)
      ecc = SQRT(ecc2)
* Calculate the new orbital angular momentum taking care to convert
* hn to units of Rsun^2/yr.
      jorb = (m1n*m2/(m1n+m2))*SQRT(hn2)*(yearsc/rsunkm)
* Determine the angle between the new and old orbital angular
* momentum vectors.
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      mu = ACOS(cmu)
* Calculate the components of the velocity of the new centre-of-mass.
 90   continue
      if(ecc.le.1.0)then
* Calculate the components of the velocity of the new centre-of-mass.
         mx1 = vk*m1n/(m1n+m2)
         mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
         vs(1) = mx1*ctheta*cphi + mx2*salpha
         vs(2) = mx1*stheta*cphi + mx2*calpha
         vs(3) = mx1*sphi
      else
* Calculate the relative hyperbolic velocity at infinity (simple method).
         sep = r/(ecc-1.d0)
*        cmu = SQRT(ecc-1.d0)
*        mu = ATAN(cmu)
         mu = ACOS(1.d0/ecc)
         vr2 = gmrkm*(m1n+m2)/sep
         vr = SQRT(vr2)
         vs(1) = vr*SIN(mu)
         vs(2) = vr*COS(mu)
         vs(3) = 0.d0
         ecc = MIN(ecc,99.99d0)
      endif
*
 95   continue
* Save some information about the colculation of the natal kick velocity
c      if(kw.eq.13)then
c         WRITE(33,*)vk,phi,theta,ecc,sep,mfin,mrem,m1n
c      elseif(kw.eq.14)then
c         WRITE(44,*)vk,phi,theta,ecc,sep,mfin,mrem,m1n
c      endif
*
      RETURN
      END
***
