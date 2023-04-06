***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      integer kw,testflag
      real*8 lum,r,mt,mc,rl,z,teff,alpha,lgdms
      real*8 ad11,ad12,ad13,ad14,ad15,ad21,ad22,ad23,ad24,ad25,ad26
      real*8 ad01,ad02,ad03,ad04,ad31,ad32,ad33,ad34,sigmaW
      real*8 gamma0,edfac,FOURPI,tem,ind,Xh
      real*8 dml,dms,dmt,p0,x,xx,mew,lum0,kap
      real*8 V2,ZsunM,eddington
      parameter(lum0=7.0d+04,kap=-0.5d0,V2=2.6d0,ZsunM=0.02)
      external eddington
*
*      windflag = 0 !BSE=0, startrack08=1, vink=2, vink+LBV for all
*      stars=3., ugonbody6=4
* Must be one of these values or mlwind will cause problem with code,
* i.e. mlwind not set (see last line of main if statement...).
    
      if(windflag.eq.0)then
* BSE
*
* Calculate stellar wind mass loss.
*
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD.
         dms = 0.d0
         if(lum.gt.4000.d0)then
            x = MIN(1.d0,(lum-4000.d0)/500.d0)
            dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
            alpha = 0.5d0
            dms = dms*(z/zsun)**(alpha)
         endif
         if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
            dml = neta*4.0d-13*r*lum/mt
            if(rl.gt.0.d0) dml =
     &         dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641,
* for high pulsation periods on AGB.
            if(kw.eq.5.or.kw.eq.6)then
               p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
               p0 = 10.d0**p0
               p0 = MIN(p0,2000.d0)
               dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
               dmt = 10.d0**dmt
               dmt = 1.d0*MIN(dmt,1.36d-09*lum)
               dml = MAX(dml,dmt)
            endif
            if(kw.gt.6)then
               dms = MAX(dml,1.0d-13*hewind*lum**(3.d0/2.d0))
            else
               dms = MAX(dml,dms)
               mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* reduced WR-like mass loss for small H-envelope mass
               if(mew.lt.1.d0)then
                  dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
                  dms = MAX(dml,dms)
               endif
* LBV-like mass loss beyond the Humphreys-Davidson limit.
               x = 1.0d-5*r*sqrt(lum)
               if(lum.gt.6.0d+05.and.x.gt.1.d0)then
                  dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
                  dms = dms + dml
               endif
            endif
         endif
*
         mlwind = dms
      elseif(windflag.eq.1)then
* StarTrack (Beclzynski+08)
*
* Calculate stellar wind mass loss.
*
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD, with no luminosity limit
* according to Belczynsk+08 pp. 174.
*
* This may not be what is actually assumed in StarTrack (see the windf1 function).
*
*
         dms = 0.d0
         if(lum.gt.4000.d0.or.(kw.ge.0.and.kw.le.1))then
            if(lum.gt.4000.d0)then
               x = MIN(1.d0,(lum-4000.d0)/500.d0)
            else
               x = 0.1d0/500.d0
            endif !or is it simply x = Min(1, lum/500)?
            dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
            alpha = 0.5d0
            dms = dms*(z/zsun)**(alpha)
         endif
         if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
            dml = neta*4.0d-13*r*lum/mt
            if(rl.gt.0.d0) dml =
     &         dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641,
* for high pulsation periods on AGB.
            if(kw.eq.5.or.kw.eq.6)then
               p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
               p0 = 10.d0**p0
               p0 = MIN(p0,2000.d0)
               dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
               dmt = 10.d0**dmt
               dmt = 1.d0*MIN(dmt,1.36d-09*lum)
               dml = MAX(dml,dmt)
            endif
            if(kw.gt.6)then
               dms = MAX(dml,1.0d-13*lum**(3.d0/2.d0)) !hewind here for KH06, not included for StarTrack...
            else
               dms = MAX(dml,dms)
               mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* reduced WR-like mass loss for small H-envelope mass
               if(mew.lt.1.d0)then
                  dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
                  dms = MAX(dml,dms)
               endif
* LBV-like mass loss beyond the Humphreys-Davidson limit.
               x = 1.0d-5*r*sqrt(lum)
               if(lum.gt.6.0d+05.and.x.gt.1.d0)then
                  dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
                  dms = dms + dml
               endif
            endif
         endif
*
         mlwind = dms
      elseif(windflag.ge.2.and.windflag.le.3) then
* Vink winds etc according to as implemented following
* Belczynski, Bulik, Fryer, Ruiter, Valsecchi, Vink & Hurley 2010.
*
* Firstly implement BSE 'old' winds that cover all other stars not
* accounted for by Vink winds (see Belczynski+09). Then implement
* Vink et al. winds.
*
* We also include the option for a variable metallicity-dependent mass
* loss parameter which eddlimflag is set, which makes the metallicity
* dependence become weaker as the star approaches the electron-scattering
* Eddington limit (Grafener & Hamann 2008, Giacobbo et al. 2018)
*
         teff = 1000.d0*((1130.d0*lum/(r**2.d0))**(1.d0/4.d0))
         dms = 0.d0
         if(lum.gt.4000.d0)then
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD after OB stars accounted for.
            x = MIN(1.d0,(lum-4000.d0)/500.d0)
            dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
            alpha = 0.5d0
            dms = dms*(z/zsun)**(alpha)
            testflag = 1
         endif
         if(kw.ge.2.and.kw.le.6)then
* 'Reimers' mass loss
            dml = neta*4.0d-13*r*lum/mt
            if(rl.gt.0.d0) dml =
     &         dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641,
* for high pulsation periods on AGB.
            if(kw.eq.5.or.kw.eq.6)then
               p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
               p0 = 10.d0**p0
               p0 = MIN(p0,2000.d0)
               dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
               dmt = 10.d0**dmt
               dmt = 1.d0*MIN(dmt,1.36d-09*lum)
               dml = MAX(dml,dmt)
            endif
            dms = MAX(dms,dml)
         endif
* Apply Vink, de Koter & Lamers (2001) OB star winds.
* Next check if hot massive H-rich O/B star in appropriate temperature ranges.
         if(teff.ge.12500.and.teff.le.25000)then
            if(eddlimflag.eq.0) alpha = 0.85d0
            dms = -6.688d0 + 2.210d0*LOG10(lum/1.0d+05) -
     &            1.339d0*LOG10(mt/30.d0) - 1.601d0*LOG10(1.3d0/2.d0) +
     &            alpha*LOG10(z/zsun) + 1.07d0*LOG10(teff/2.0d+04)
            dms = 10.d0**dms
            testflag = 2
         elseif(teff.gt.25000.)then
*        Although Vink et al. formulae  are only defined until Teff=50000K,
*        we follow the Dutch prescription of MESA, and extend to higher Teff
             dms = -6.697d0 + 2.194d0*LOG10(lum/1.0d+05) -
     &            1.313d0*LOG10(mt/30.d0) - 1.226d0*LOG10(2.6d0/2.d0) +
     &            alpha*LOG10(z/zsun) +0.933d0*LOG10(teff/4.0d+04) -
     &            10.92d0*(LOG10(teff/4.0d+04)**2)
       dms = 10.d0**dms
       testflag = 2
         endif

         if((windflag.eq.3.or.kw.ge.2).and.kw.le.6)then
* LBV-like mass loss beyond the Humphreys-Davidson limit.
* Optional flag (windflag=3) to use for every non-degenerate star
* past the limit, rather than just for giant, evolved stars
            x = 1.0d-5*r*sqrt(lum)
            if(lum.gt.6.0d+05.and.x.gt.1.d0)then
               if(eddlimflag.eq.0) alpha = 0.d0
               dms = 1.5d0*1.0d-04*((z/zsun)**alpha)
               testflag = 3
            endif
         elseif(kw.ge.7.and.kw.le.9)then !WR (naked helium stars)
* If naked helium use Hamann & Koesterke (1998) WR winds reduced by factor of
* 10 (Yoon & Langer 2005), with Vink & de Koter (2005) metallicity dependence
            if(eddlimflag.eq.0) alpha = 0.86d0
            dms = 1.0d-13*(lum**1.5d0)*((z/zsun)**alpha)
            testflag = 4
         endif
*
         mlwind = dms
         
      elseif(windflag.eq.4) then
* Calculate stellar wind mass loss.
***
* Calculate the effective temperature (tem) by the simple formula:
*     L = 4*pi*T^4*sigmaW
             sigmaW = 5.67d0*10.d0**(-5.d0)*(6.96d0*10.d0**10.d0)**2.d0
     &        /(3.84d0*10.d0**33.d0)
               FOURPI = 2.d0*ACOS(-1.d0)
              tem = (lum/(FOURPI*r**2.d0*sigmaW))**(1.d0/4.d0)
*
* Exponent of the dependence on metallicity: (Z/ZsunM)^ind
* Dipendence on eddingtin factor taken from Chen et al., MNRAS,
* 452, 2015.
* Eddington factor is given by eddinton()
              ind = 0.85d0
              edfac = eddington(mt,lum,kw)
              if(2.d0/3.d0.lt.edfac.and.edfac.lt.1.d0)then
                 ind = (2.45d0 - 2.4d0*edfac)
              elseif(edfac.ge.1.d0)then
                 ind = 0.05d0
              endif
*
***
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD.
*
              dms = 0.d0
              if(lum.gt.4000.d0)then
                 x = MIN(1.d0,(lum-4000.d0)/500.d0)
                 dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
                 dms = dms*(z/ZsunM)**(ind)
              endif
*
* Stellar winds for O/B stars.
*
              if(kw.ge.0 .and. kw.le.1)then
*
* Eqs 25 and 24 in Vink et al. 2001, A&A, 369, 574.
*
                  if(12500.d0.le.tem .and. tem.le.25000.d0)then
                     ad11 = 2.21d0*Log10(lum/(10.d0**5.d0))
                     ad12 = -1.339d0*Log10(mt/30.d0)
* V is the ration of wind velocity at inf to escape velosity.
                     ad13 = -1.601d0*Log10(1.3d0/2.d0)
                     ad14 = ind*Log10(z/ZsunM)
                     ad15 = 1.071*Log10(tem/20000.d0)
* logarithm of the mass loss rate.
                     lgdms = -6.688d0 + ad11 + ad12 + ad13+ad14+ad15
                     dms = 10.d0**lgdms
                  elseif(25000.d0.lt.tem .and. tem.le.50000.d0)then
                     ad21 = 2.194d0*Log10(lum/(10.d0**5.d0))
                     ad22 = -1.313d0*Log10(mt/30.d0)
* V is the ration of wind velocity at inf to escape velosity.
                     ad23 = -1.226d0*Log10(V2/2.d0)
                     ad24 = ind*Log10(z/ZsunM)
                     ad25 = 0.933d0*Log10(tem/40000.d0)
                     ad26 = -10.92d0*(Log10(tem/40000.d0))**2.d0
* logarithm of the mass loss rate
                     lgdms = -6.697d0 + ad21 + ad22+ad23+ad24+ad25+ad26
                     dms = 10.d0**lgdms
                  endif
              endif
*
***
*
              if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
                 dml = neta*4.0d-13*r*lum/mt
*
* Check for any tidally enhanced mass loss in binary systems (optional):
* see Tout & Eggleton, MNRAS, 1988, 231, 823.
*
                 if(rl.gt.0.d0) dml=dml*(1.d0+bwind*
     &          (MIN(0.5d0,(r/rl)))**6)

                 dms = max(dms,dml)
*
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641,
* for high pulsation periods on AGB.
*
                 if(kw.eq.6)then
                    p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
                    p0 = 10.d0**p0
                    p0 = MIN(p0,2000.d0)
                    dmt = -11.4d0+0.0125d0*(p0-100.d0*
     &              MAX(mt-2.5d0,0.d0))
                    dmt = 10.d0**dmt
                    dmt = 1.d0*MIN(dmt,1.36d-09*lum)
                    dms = MAX(dms,dmt)
                    endif
                 if(kw.gt.6)then
*
***
* WR-like mass loss from Belczynski 2010 plus the dependence on 
* metallicity (Giacobbo et al. 2017).
*
                    dms = (10.d0**(-13.d0))*(lum**1.5d0)*(z/ZsunM)**ind
*
                 endif
*
                 xx = 1.0d-5*r*sqrt(lum)
                 if(lum.gt.6.0d+05.and.xx.gt.1.d0)then
*
* LBV-like mass loss beyond the Humphreys-Davidson limit
* and beyond the MS (it depends on Z)
*
                    dms = (z/ZsunM)**ind*1.5d0*10.d0**(-4.d0)
*
                 endif
              endif
*
              mlwind = dms
      endif

      return
      end
***
