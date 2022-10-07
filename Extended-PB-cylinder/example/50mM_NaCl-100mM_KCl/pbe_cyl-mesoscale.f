      program pbe

      implicit none

c********************************************************************c
c  numerical solution of Poisson-Boltzmann equation for ions at the  c
c  water/hydrophobic interface.
c  Surface has a Charge Density of sigma                             c
c********************************************************************c
c
c********************************************************************c
c  Modification to account for ion-competition around RNA            c
c  PB solved in cylindrical coordinates                              c
c  Solve for 3 PMFs and 3 concentration                              c
c********************************************************************c
c
c input pmf: units are kT and nm
c input: sigma in e/nm
c
c *** number of gridpoints ***
      integer ngrd
      integer ngrd2
      parameter (ngrd =30000)
      parameter (ngrd2=15000)

c *** Startvalue of Phi_ion(1)
c      real*16 SV
c      parameter(SV= -10.0)

c *** number of concentrations scanned ***
      integer n_c

c *** print level ***
      integer iprint

c *** free energies ***
c *** Add E3_pmf to include the cation competition (Sergio)***
      real*16 E1_pmf(ngrd), E2_pmf(ngrd), E3_pmf(ngrd)

c *** Add E4_pmf=cation1-nucleobase and E5_pmf=cation2 nucleobase
c *** to include the mesoscale approach (Sergio)***
      real*16 E4_pmf(ngrd), E5_pmf(ngrd), E6_pmf(ngrd), E7_pmf(ngrd)
c *** intrinsic, ionic, and total potential ***
      real*16 Phi_int(ngrd), Phi_tot(ngrd), Phi_ion(ngrd)

c *** derivatives and PB, el. field ***
      real*16 PB(ngrd) , el_field(ngrd)

c *** Density Dependant Epsilon ***
      real*16 epsilon(ngrd), r_eps(ngrd)
      real*16 PB1(ngrd), PB2(ngrd)
     
c *** distances from data files
      real*16 r_phi(ngrd), r_pmf(ngrd)

c *** concentration profiles, surface excess
      real*16 c1(ngrd), c2(ngrd), Gamma_1, Gamma_2,Gamma_3
      real*16 c3(ngrd)
      real*16 g1(ngrd), g2(ngrd), g3(ngrd)

c *** surface tension
      real*16 surf_tens, dc, surf_tens_rrn
      real*16 st_1, st_2, st_3, st_4, st_5
      real*16 st1, st2, st3, st4, st5

c *** Gouy-Chapman potential ***
      real*16 phi_gc, phi_gc0,gc_fac
      real*16 kappa_dh
      real*16 sigma, sigma_el, sigma_el0, sigma_fac
      real*16 sigma_tot

c *** gibbs dividing surface, depletion layer thickness ***
      integer r_gds
      real*16 r_depl

c *** radius
      real*16 radius

c *** distance between the plates (in nm) ***
      real*16 d, dd, dp
      integer n_middle

c *** Zeta parameter Langmuir 2013, 29, 2602−2614 eq. (3)
c      real*16 zeta
c ******* define units **********
c *** c:                mol/l ***
c *** r:                nm    ***
c *** Phi:              V     ***
c *** q:                e     ***
c *** E:                kT    ***
c *** rho_el:         e*mol/l ***
c *** sigma:          e/nm^2  ***
c *** sigma_el:       e/nm^2  ***
c *** EF:               V/nm  ***
c *** T:                300 K ***
c********************************
c *** phi_fac converts phi into volt ***
c *** q_fac = e_0/kT (300K) ***
c *** phi_fac = (l*eps_0*eps_w)/(e*N_A*nm^2) *** 
c *** parameter (k_B = 1.380662e-23 J/K) ***
c *** parameter (e = 1.6021892e-19 C) ***
c *** parameter (eps_0 = 8.8541878e-12 C^2/Jm) ***
c *** parameter (eps_w = 80 ) ***
c *** parameter (N_A = 6.022045e23 1/mol) ***
c *** parameter (gamma_fac = 1000 nm*k*N_A*T/l) (mN/m^2) ***
c *** parameter (sigma_fac = N_A/1e24) ***
c *** parameter (gc_fac = e/(eps0* eps_w * 10^-9 ) ***
c *** parameter (pol_fac = 4 pi eps_0 V^2/nm^2 Ang^3 / k_B*T) ***
c *** parameter (number_fac = 2 pi N_A/10^24 = 3.782477) ***

      real*16 q_fac, phi_fac, gamma_fac, pol_fac
      real*16 number_fac
      parameter (phi_fac = 0.1397)
      parameter (q_fac = 38.681666)
      parameter (gamma_fac = 2.494e0)
      parameter (sigma_fac = 0.6022)
      parameter (gc_fac = 2.2619e-1)
      parameter (pol_fac = 0.02686)
      parameter (number_fac = 3.78247755)
      

c *** distance element ***
      real*16 dr
      parameter (dr=0.001)

c *** bulk concentration, charge , polarizabilities ***
      real*16 c1_0, c2_0, q1, q2, c1_00, c1_01, pol1, pol2

c *** bulk concentration, and charge cation type 2, and aux***
      real*16 c3_0, q3, c3_01

c *** iteration ***
      real*16 fail_low, fail_high, eps_conv
      integer failstep, failmode, i_shot, n_shot
      integer n_conv, i_conv
      integer index
      parameter (n_shot=200)

c *** io units **
      integer outfile, outfile2, outfile3, pmf1file, pmf2file, phifile
      integer outfile4, outfile5, outfile6, outfile7
      integer outfile10, outfile11, outfile12, outfile13, outfile14
      integer outfile15, outfile16, outfile17, outfile18, pmf3file
      integer pmf4file, pmf5file,pmf6file, pmf7file
      integer parmfile, epsfile
      parameter (outfile=10, outfile2=11) 
c *** Include pmf3file=28, and above define the pmf3file
      parameter (pmf1file=12, pmf2file=13, phifile=14)   
      parameter (parmfile=15)
      parameter (outfile3=16)
      parameter (outfile4=17)
      parameter (outfile5=18)
      parameter (outfile6=19)
      parameter (outfile11=20, outfile12=21, outfile13=22)
      parameter (outfile14=23, outfile15=24, outfile16=25)
      parameter (outfile17=26, outfile18=27, pmf3file=28)
      parameter (pmf4file=29, pmf5file=30, pmf6file=31, pmf7file=32)

c *** etc ***
      integer i, j, k
      real*16 deltaE1, deltaE2


c *** Numbers to calculate the energies ***
      real*16 deltaE3
c *** Numbers to calculate the energies mesoscale***
      real*16 deltaE4, deltaE5, deltaE6, deltaE7, zeta


c *** stuff for the restricted densities 
c *** cn1, cn2, cn3, nen1, nen2, nen, dia
      real*16 cn1(ngrd), cn2(ngrd), cn3(ngrd)
      real*16 nen(ngrd), nen1(ngrd), nen2(ngrd)
      real*16 dia1, dia2, dia3

      open (outfile,file='pbe.out',status='replace')
      open (outfile2,file='pbe.out2',status='replace')
      open (outfile3,file='TDensity0_01M.dat',status='replace')
      open (outfile4, file='TDensity0_1M.dat',status='replace')
      open (outfile5,file='TDensity1M.dat',status='replace')
      open (outfile10, file='Density0_1M.dat',status='replace')
      open (outfile11, file='Density0_2M.dat',status='replace')
      open (outfile12, file='Density0_3M.dat',status='replace')
      open (outfile13, file='Density0_4M.dat',status='replace')
      open (outfile14, file='Density0_5M.dat',status='replace')
      open (outfile15, file='Density0_6M.dat',status='replace')
      open (outfile16, file='Density0_7M.dat',status='replace')
      open (outfile17, file='Density0_8M.dat',status='replace')
      open (outfile18, file='Density0_9M.dat',status='replace')
      open (parmfile,file='pbe.parm',status='old')
      open (epsfile,file='Epsilon_z.dat',status='old')
      open (outfile6,file='GC.dat',status='replace')


      read (parmfile,*)
      read (parmfile,*) q1, q2, q3
      read (parmfile,*)
      read (parmfile,*) c1_00, c1_01, c3_01
      read (parmfile,*)
      read (parmfile,*) n_c
      read (parmfile,*)
      read (parmfile,*) iprint
      read (parmfile,*)
      read (parmfile,*) r_gds, r_depl
      read (parmfile,*)
      read (parmfile,*) sigma
      read (parmfile,*)
      read (parmfile,*) radius
      read (parmfile,*)
      read (parmfile,*) d
      read (parmfile,*)
      read (parmfile,*) dia1, dia2, dia3
      read (parmfile,*)
      read (parmfile,*) zeta

      if (iprint.ge.1) then
        write (outfile,'(a60, a60)') 
     1 "#     c1_0          Phi(1) Phi(2)   Gamma_1     Gamma_2     ",
     2 "  surf_tens     surf_tens_rrn                               "
      end if

      if (iprint.ge.1) then
        write (outfile2,'(a60, a60)')
     1 "#     dr*i          c1           c2          c3   g1(i)  ",
     2 "      g2(i)         g(3)         Gamma_1     Gamma_2     ",
     3 "#  Gamma_3 "
      end if

      if (iprint.ge.1) then
        write (outfile3,'(a60, a60)')
     1 "#     dr*i          c1           c2     Phi(1)  "
      end if

      if (iprint.ge.1) then
        write (outfile4,'(a60, a60)')
     1 "#     dr*i          c1           c2      Phi(1)   "
      end if

      if (iprint.ge.1) then
        write (outfile5,'(a60, a60)')
     1 "#     dr*i          c1           c2    Phi(1)     "
      end if

      if (iprint.ge.1) then
         write (outfile6,'(a60, a60)')
     1 "sigma_el0, sigma_el, kappa_dh, Phi_gc0, Phi_gc"
      end if

      dc = (c1_01 - c1_00) / float(n_c) 
      write (*,*) 'q1 = ', q1, '   q2 = ', q2
      write (*,*) 'c1:  ', c1_00, ' to ', c1_01, n_c, '  steps'


c *** read potential ***
      open (phifile,file='potential.dat',status='old')
      do i = 1, ngrd
        read (phifile,*) r_phi(i), Phi_int(i)
      enddo

c *** read pmf ***
      open (pmf1file,file='E1_pmf.dat',status='old')
      do i = 1, ngrd
        read (pmf1file,*) r_pmf(i), E1_pmf(i)
      enddo

c *** read pmf ***
      open (pmf2file,file='E2_pmf.dat',status='old')
      do i = 1, ngrd
        read (pmf2file,*) r_pmf(i), E2_pmf(i)
      enddo


c *** read pmf3 - competitive-backbone***
      open (pmf3file,file='E3_pmf.dat',status='old')
      do i = 1, ngrd
        read (pmf3file,*) r_pmf(i), E3_pmf(i)
      enddo

c *** read pmf4 - buffer-nucleobase ***
      open (pmf4file,file='E4_pmf.dat',status='old')
      do i = 1, ngrd
        read (pmf4file,*) r_pmf(i), E4_pmf(i)
      enddo

c *** read pmf5 - competitive-nucleobase ***
      open (pmf5file,file='E5_pmf.dat',status='old')
      do i = 1, ngrd
        read (pmf5file,*) r_pmf(i), E5_pmf(i)
      enddo

c *** read pmf6 - competitive-nucleobase O6 ***
      open (pmf6file,file='E6_pmf.dat',status='old')
      do i = 1, ngrd
        read (pmf6file,*) r_pmf(i), E6_pmf(i)
      enddo

c *** read pmf7 - competitive-nucleobase ***
      open (pmf7file,file='E7_pmf.dat',status='old')
      do i = 1, ngrd
        read (pmf7file,*) r_pmf(i), E7_pmf(i)
      enddo

c *** read Epsilon(z) ***
      do i = 1, ngrd
        read (epsfile,*) r_eps(i), epsilon(i)
      enddo

c      do i=1, ngrd
c         write(*,*) epsilon(i), i
c      enddo

c *** correct long tail of pmfs
      do i = 3001, ngrd
        E1_pmf(i) = 0.0
	E2_pmf(i) = 0.0
      enddo

c *** initialize surface tension ***
      surf_tens = 0.0d0


c *** loop over concentrations ****
c -@Sergio: Here initiate the large loop!!!
c      do  k = 1, n_c
        
c        c1_0 = c1_00 + k*((c1_01-c1_00)/n_c)
c        c2_0 = - c1_0 * q1/q2

c *** @Sergio Modification to include competition of ions
c     C2 is fixed, C3 varies for competition, and C1 is set
c     for a neutral system
      do  k = 1, n_c
        c2_0 = c1_01
        c3_0 = c1_00 + k*((c3_01-c1_00)/n_c)
        c1_0 = - (c2_0 * q2/q1 + c3_0 * q3/q1)


        if ((q1*c1_0+q2*c2_0 + q3*c3_0)**2.ge.1e-9) then
           write (*,'(a20, f12.8, f12.8)') 'System is  charged!!!',
     1      q1*c1_0,q2*c2_0, q3*c3_0
        end if      

c *** initialize ***
        do i = 1, ngrd
          Phi_tot(i) = 0.0
          Phi_ion(i) = 0.0
          PB(i) = 0.0
          PB1(i) = 0.0
          PB2(i) = 0.0
          el_field(i) = 0.0
          c1(i) = c1_0
          c2(i) = c2_0
        enddo


c *** shooting from some phi(0), we want to get asymptotically ***  
c *** to phi (inf) = 0 ***
c *** assume phi(0). start with phi(0) = -10 (surely too high) ***
c *** shot failed if |phi| > 10 ***
c *** or if phi changes sign for z > 15 nm ***
c *** E(1) is set by sigma ***
c *** Calculte Phi(2) from Phi(1) + dr * E(1) ***
c *** Use eps_w =80 for BC ***

        i_conv=0
        n_conv=5000
        eps_conv = 1.q-8
        Phi_ion(1) = -10.0
        el_field(1) = - sigma * gc_fac

        PB1(1) = -sigma * gc_fac*(radius*dr) 
        PB2(1) = -sigma  * gc_fac 

        Phi_ion(2) = Phi_ion(1) + el_field(1)*dr

        fail_low = Phi_ion(1)
        fail_high = -fail_low
        failmode = 0
        i_shot = 0
        index =0

c *** continue after shot has failed ***
 777    continue
        i_shot = i_shot + 1

c *** update intervall ***
        if (failmode .eq. 1) fail_low = Phi_ion(1)
        if (failmode .eq. 2) fail_high = Phi_ion(1)

c *** update Psi at the wall ***
        Phi_ion(1) = fail_high + 0.5*(fail_low-fail_high)
        Phi_ion(2) = Phi_ion(1) + el_field(1)*dr
        if (failmode .eq. 0) then
           Phi_ion(1) = -10
           Phi_ion(2) = -10
        end if


c *** reset the failmode ***
        failmode = 0

        do i = 3, ngrd

          deltaE1 = (E1_pmf(i-1) - E1_pmf(ngrd)) + q1*q_fac*Phi_ion(i-1) 
          cn1(i-1) = c1_0 * exp(-deltaE1)

c ***  Include the mesoscale approach for the   
        deltaE2 = (E2_pmf(i-1) - E2_pmf(ngrd)) + q2*q_fac*Phi_ion(i-1)
        deltaE4 = (E4_pmf(i-1)- E4_pmf(ngrd)) + q2*q_fac*Phi_ion(i-1) 
c ***  Include the mesoscale - O6 - site
        deltaE6 = (E6_pmf(i-1)- E6_pmf(ngrd)) + q2*q_fac*Phi_ion(i-1) 
c ***     cn2(i-1) = c2_0 * exp(-deltaE2): 2 backbone - 1 N7 - 1 O6 (0.5 - 0.25 - 0.25)
          cn2(i-1) = c2_0 *zeta* exp(-deltaE2) + 
     1            c2_0 *(0.25)* exp(-deltaE4)+
     1            c2_0 *(0.25)* exp(-deltaE6)

          deltaE3 = (E3_pmf(i-1) - E3_pmf(ngrd)) + q3*q_fac*Phi_ion(i-1)
          deltaE5 = (E5_pmf(i-1)- E5_pmf(ngrd)) + q3*q_fac*Phi_ion(i-1)
c ***  Include the mesoscale - O6 - site 
        deltaE7 = (E7_pmf(i-1)- E7_pmf(ngrd)) + q2*q_fac*Phi_ion(i-1)  
c ***     cn3(i-1) = c3_0 * exp(-deltaE3): 2 backbone - 1 N7 - 1 O6 (0.5 - 0.25 - 0.25)
          cn3(i-1) = c3_0 *zeta* exp(-deltaE3) + 
     1            c3_0 *(0.25)* exp(-deltaE5)+
     1            c3_0 *(0.25)* exp(-deltaE7)

       

c     Modification for finite concentration - Sergio Cruz
c     Neue Dichten: cn1, cn2, cn3, nen1, nen2, nen, dia                                                                                                       
          nen1(i-1) = sqrt(2.0) + dia2*dia2*dia2*(cn2(i-1)-c2_0)
          nen2(i-1) = dia1*dia1*dia1*(cn1(i-1)-c1_0) + 
     1      dia3*dia3*dia3*(cn3(i-1)-c3_0)
          nen(i-1) = nen1(i-1) + nen2(i-1)

          c1(i-1) = (cn1(i-1)*sqrt(2.0))/(nen(i-1))

          c2(i-1) = (cn2(i-1)*sqrt(2.0))/(nen(i-1))
          c3(i-1) = (cn3(i-1)*sqrt(2.0))/(nen(i-1))
c     Finish of the modification- Sergio Cruz
 
          PB(i-1) =-(radius*dr+(i-1)*dr)*phi_fac*(q1*c1(i-1)
     1      +q2*c2(i-1)+q3*c3(i-1))

          PB1(i-1)=  PB1(i-2) + PB(i-1) * dr

          el_field(i-1) = PB1(i-1)/(radius*dr+(i-1)*dr) 

          PB2(i-1) = PB1(i-1)/(radius*dr+(i-1)*dr)
          
          Phi_ion(i) = Phi_ion(i-1) +  dr * PB2(i-1)

    
c *** check for failure ***
c *** failmode:   
c     1  phi too low
c     2  phi too high

          if (Phi_ion(i) .le. -50 ) then
            failstep = i
            failmode = 1
            index = 1
             if (i_shot.le.n_shot) goto 777
          end if
          
          if (Phi_ion(i) .ge. 50 ) then
            failstep = i
            failmode = 2
            index =2
             if (i_shot.le.n_shot) goto 777
          end if


c *** Change of sign ***
c *** Only usefull is PMF = 0 ***
c *** Here: For z > 15 nm ***

          if(i .ge. 15000) then
             if (Phi_ion(i)*Phi_ion(i-1) .lt .0.0) then
                if (i_shot .eq. 1) then
                   write(*,*) "WARNING: Step1 phi is changing sign!"
                   write(*,*) "Change start value of potential!"
                end if
                if (Phi_ion(1).le.0) then
                   failstep = i
                   failmode = 2
                   index =3
                   if (i_shot.le.n_shot) goto 777
                else
                   failstep = i
                   failmode = 1
                   index =4
                   if (i_shot.le.n_shot) goto 777
                end if
                
             end if
          end if
       enddo 


c *** 2nd check for failure ***
          if (Phi_ion(ngrd2) .le. -eps_conv ) then
            failstep = i
            failmode = 1
            write(*,*) "phi too high without breaking bounds"
            write(*,*) "at the end of the 0.5 grid"
            write(*,*) i_shot
            write(*,'(a18,f12.6)') "wall potential: ",Phi_ion(ngrd)
             if (i_shot.le.n_shot) goto 777
          end if

          if (Phi_ion(ngrd2) .ge. eps_conv ) then
            failstep = i
            failmode = 2
            write(*,*) "phi too low without breaking bounds"
            write(*,*) "at the end of the 0.5 grid"
            write(*,*) i_shot
            write(*,'(a18,f12.6)') "wall potential: ",Phi_ion(ngrd)
             if (i_shot.le.n_shot) goto 777
          end if

c *** 3rd check for failure ***
          if (Phi_ion(ngrd) .le. -eps_conv ) then
            failstep = i
            failmode = 1
            write(*,*) "phi too high without breaking bounds"
            write(*,*) "at the end of the grid"
            write(*,*) i_shot
            write(*,'(a18,f12.6)') "wall potential: ",Phi_ion(ngrd)
             if (i_shot.le.n_shot) goto 777
          end if

          if (Phi_ion(ngrd) .ge. eps_conv ) then
            failstep = i
            failmode = 2
            write(*,*) "phi too low without breaking bounds"
            write(*,*) "at the end of the grid"
            write(*,*) i_shot
            write(*,'(a18,f12.6)') "wall potential: ",Phi_ion(ngrd)
             if (i_shot.le.n_shot) goto 777
          end if


c ***        c1(1) = 0.0 
c ***        c2(1) = 0.0

        deltaE1 = (E1_pmf(1) - E1_pmf(ngrd)) + q1*q_fac*Phi_ion(1)                                                                                   
c ****  c1(1) = c1_0 * exp(-deltaE1)

        deltaE3 = (E3_pmf(i-1) - E3_pmf(ngrd)) + q3*q_fac*Phi_ion(i-1) 
        c3(i-1) = c3_0* exp(-deltaE3) 


        deltaE2 = (E2_pmf(1) - E2_pmf(ngrd)) + q2*q_fac*Phi_ion(1)
        c2(1) = c2_0 * exp(-deltaE2)

c *** debug output ***
        if (iprint.ge.3) then
          write (*,*) '2converged w/o pol'
        end if

c *** *********************************************************************************************


c *** calculate the Gouy-Chapman potential ***

      sigma_el0 = 0.0
      sigma_el = 0.0
      sigma_tot = 0.0
      do i = 1, ngrd

       sigma_el0 = sigma_el0 + 
     1            dr * q1 * c1_0 * exp(-(E1_pmf(i) - E1_pmf(ngrd)))
       sigma_el0 = sigma_el0 +
     1            dr * q2 * c2_0 * exp(-(E2_pmf(i) - E2_pmf(ngrd)))

       sigma_el = sigma_el +
     1            dr * q1 * c1_0 * exp(-(E1_pmf(i) - E1_pmf(ngrd) +
     2            q1*q_fac*Phi_ion(i)))
       sigma_el = sigma_el +
     1            dr * q2 * c2_0 * exp(-(E2_pmf(i) - E2_pmf(ngrd) + 
     2            q1*q_fac*Phi_ion(i)))

      enddo

      sigma_el0 = sigma_el0 * sigma_fac
      sigma_el = sigma_el * sigma_fac
      sigma_tot = sigma_el + sigma
      kappa_dh = 3.246 * sqrt(c1_0)
      Phi_gc0 = sigma_el0 * gc_fac / kappa_dh
      Phi_gc = sigma_el * gc_fac / kappa_dh


c *** converged, calculate surface excess ***
c *** avoid the tail ***
c *** r_gds is the pos. of the  gibbs dividing surface ***
c *** check ... ***

      Gamma_1 = 0.0
      Gamma_2 = 0.0
      Gamma_3 = 0.0


c *** count the density outside the GDS ****
c      do i = 10, r_gds - 1
c
c        Gamma_1 = Gamma_1 + dr * c1(i) 
c        Gamma_2 = Gamma_2 + dr * c2(i) 
c        g1(i) = Gamma_1
c        g2(i) = Gamma_2

c      enddo

      do i = r_gds, 15000

        Gamma_1 = Gamma_1+
     1            number_fac*dr* (c1(i) - c1_0)*(radius*dr+(i-1)*dr)
        Gamma_2 = Gamma_2+
     1            number_fac*dr* (c2(i) - c2_0)*(radius*dr+(i-1)*dr)
        Gamma_3 = Gamma_3+
     1            number_fac*dr* (c3(i) - c3_0)*(radius*dr+(i-1)*dr)
        g1(i) = Gamma_1
        g2(i) = Gamma_2
        g3(i) = Gamma_3

      enddo

c *** surface tension ***

       surf_tens = surf_tens - dc * (Gamma_1 + Gamma_2) / c1_0
     1            * gamma_fac

cc *** use this part for inclusion of DH bulk ion activity ***
c      surf_tens = surf_tens - dc * (Gamma_1 + Gamma_2) / c1_0
c     1            * gamma_fac +
c     2         dc * (Gamma_1 + Gamma_2) * 0.586 * sqrt(c1_0)
c ************************************************************

c *** surface tension following roland ***

      surf_tens_rrn = 0.0
      st1 = 0
      st2 = 0
      st3 = 0
      st4 = 0
      st5 = 0

      do i = 1, r_gds-1

        st_1 = (q1*c1(i) + q2*c1(i)) * q_fac * 0.5 * Phi_ion(i) 
        st_2 = c1(i) * log ( c1(i)/ c1_0 ) - (c1(i) - c1_0)
        st_3 = c2(i) * log ( c2(i)/ c1_0 ) - (c2(i) - c1_0)           
        st_4 = (E1_pmf(i) - E1_pmf(ngrd)) * c1(i)
        st_5 = (E2_pmf(i) - E2_pmf(ngrd)) * c2(i)

        st1 = st1 + st_1
        st2 = st2 + st_2
        st3 = st3 + st_3
        st4 = st4 + st_4
        st5 = st5 + st_5

        surf_tens_rrn = surf_tens_rrn +
     1                  (st_1 + st_2 + st_3 + st_4 + st_5) * dr
     2                * gamma_fac

      enddo

      do i = r_gds, ngrd/4

        st_1 = (q1*c1(i) + q2*c1(i)) * q_fac * 0.5 * Phi_ion(i) 
        st_2 = c1(i) * log ( c1(i)/ c1_0 ) - (c1(i) - c1_0)
        st_3 = c2(i) * log ( c2(i)/ c1_0 ) - (c2(i) - c1_0)           
        st_4 = (E1_pmf(i) - E1_pmf(ngrd)) * c1(i)
        st_5 = (E2_pmf(i) - E2_pmf(ngrd)) * c2(i)

        st1 = st1 + st_1
        st2 = st2 + st_2
        st3 = st3 + st_3
        st4 = st4 + st_4
        st5 = st5 + st_5

        surf_tens_rrn = surf_tens_rrn + 
     1                  (st_1 + st_2 + st_3 + st_4 + st_5) * dr
     2                * gamma_fac  

      enddo


c *** final output ***

c      if ((iprint.ge.1) .and. (k.eq.n_c)) then
c           do i = 1, ngrd
c               write (outfile2,'(12(2x,f12.6))')
c     1              dr*i, c1(i)/c1_0,c2(i)/c2_0, c3(i)/c3_0,
c     1              Phi_ion(i),
c     2              g1(i),  g2(i),  Gamma_1, Gamma_2
c           enddo
c      end if

      if ((iprint.ge.1) .and. (k.eq.n_c)) then
           do i = 1, ngrd
               write (outfile2,'(12(2x,f12.6))')
     1              dr*i, c1(i), c2(i), c3(i),
     1              Phi_ion(i), g1(i),  g2(i), g3(i),
     2              Gamma_1, Gamma_2, Gamma_3
           enddo
      end if


c *** this is the output of all concentrations... ***

      if ((iprint.ge.1) .and. (k.eq.1)) then
           do i = 1, ngrd
               write (outfile3,'(3(2x,f12.6), 1(2x,g16.8))')
     1         dr*i, c1(i)/c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if


      if ((iprint.ge.1) .and. (k.eq.10)) then
           do i = 1, ngrd
               write (outfile4,'(3(2x,f12.6), 1(2x,g16.8))')
     1              dr*i, c1(i)/c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if

      if ((iprint.ge.1) .and. (k.eq.100)) then
           do i = 1, ngrd
               write (outfile5,'(3(2x,f12.6), 1(2x,g16.8))')
     1              dr*i, c1(i)/c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if

      if ((iprint.ge.1) .and. (k.eq.20)) then
           do i = 1, ngrd
               write (outfile11,'(3(2x,f12.6), 1(2x,g16.8))')
     1              dr*i,c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if

      if ((iprint.ge.1) .and. (k.eq.30)) then
           do i = 1, ngrd
               write (outfile12,'(3(2x,f12.6), 1(2x,g16.8))')
     1              dr*i, c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if

      if ((iprint.ge.1) .and. (k.eq.40)) then
           do i = 1, ngrd
               write (outfile13,'(3(2x,f12.6), 1(2x,g16.8))')
     1              dr*i, c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if

      if ((iprint.ge.1) .and. (k.eq.50)) then
           do i = 1, ngrd
               write (outfile14,'(3(2x,f12.6), 1(2x,g16.8))')
     1              dr*i, c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if

      if ((iprint.ge.1) .and. (k.eq.60)) then
           do i = 1, ngrd
               write (outfile15,'(3(2x,f12.6), 1(2x,g16.8))')
     1              dr*i, c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if

      if ((iprint.ge.1) .and. (k.eq.70)) then
           do i = 1, ngrd
               write (outfile16,'(3(2x,f12.6), 1(2x,g16.8))')
     1              dr*i, c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if

      if ((iprint.ge.1) .and. (k.eq.80)) then
           do i = 1, ngrd
               write (outfile17,'(3(2x,f12.6), 1(2x,g16.8))')
     1              dr*i, c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if

      if ((iprint.ge.1) .and. (k.eq.90)) then
           do i = 1, ngrd
               write (outfile18,'(3(2x,f12.6), 1(2x,g16.8))')
     1              dr*i, c1_0,c2(i)/c2_0,Phi_ion(i)
           enddo
      end if


c *** End of Concentration output ***

        if (iprint.ge.1) then
          write (outfile,('(1(2x,f12.6),5(2x,g16.8))'))
     1           c1_0, Phi_ion(1),
     1           Phi_ion(2),  Gamma_1, Gamma_2, surf_tens
          write (*,*)  k, c1_0
        end if

c *** Gouy Chapman Results ***

         write (outfile6,('(7(2x,f12.6))'))
     1           c1_0, sigma_el0, sigma_el, kappa_dh, Phi_gc0,
     1           Phi_gc, sigma_tot


      enddo

      end

