!  *********************************************************
!   molecular dynamics simulation on protein binding and folding
!                                     
!   author: Hueseyin Kaya            revised by Zhirong Liu
! 
!                                           by Yongqi Huang
! 


!   SchemeO (0): chain A is moves without any influence from chain B (by ignoring interchain Q)
!   SchemeA (1): chain B is fixed to avoid unfolding, but chain A is allowed to move freely
!   SchemeB (2): both chains A and B are moving
!  
! 	Periodic boundary condiction is used to make sure chain A is near chain B.
! 
!   !  Solvation model: M. S. Cheung et al, PNAS 99, 685 (2002) 
!   note that there is a little mistake for term C in the formula
!   
!   Simulate the system to see if it first reach the native 
!   or the denatural states
! 
!   restore when the data is running out
!   ouput some conformation for further usage
! 
! 
!   scale the bending, torsion, non-bonding interactions
!   Non-native HP interactions with changable strength.

!   Note that we use the appNCS file in data format of: residue_i, residue_j, contact_type (1: native contact; 2: non-native hydrophobic interaction;  else: repulsive)
!   A gap between residue number of two chains can be specified.
!   Other columns than the first three are ignored.
!   The order and number of contact intems can be random. 

!   Add the bias and histogram for Qw, which is defined as Qw=Qb-\alpha_w*\sum_native [r_AB-r_AB^(0)]
!   The order and number of contact intems can be random. 


!  References:
!  [1] Yongqi Huang and Zhirong Liu. 
!      Kinetic advantage of intrinsically disordered proteins in coupled folding-binding process: a critical assessment of the ¡°fly-casting¡± mechanism. 
!      J. Mol. Biol. 393 (5), 1143-1159 (2009). 
!  [2] Yongqi Huang and Zhirong Liu. 
!      Smoothing molecular interactions: the ¡°kinetic buffer¡± effect of intrinsically disordered proteins. 
!      Proteins: Struct. Funct. Bioinform. 78 (16), 3251-3259 (2010).
!  [3] Yongqi Huang and Zhirong Liu. 
!      Nonnative interactions in coupled folding and binding processes of intrinsically disordered proteins. 
!      PLoS ONE 5 (11), e15375/1-10 (2010).
!  [4] Huseyin Kaya and Hue Sun Chan.
!      Solvation Effects and Driving Forces for Protein Thermodynamic and Kinetic Cooperativity: How Adequate is Native-centric Topological Modeling?
!      J. Mol. Biol. 326 (3), 911-931 (2003). 
!  [5] Zhirong Liu and Hue Sun Chan. 
!      Solvation and desolvation effects in protein folding: native flexibility, kinetic cooperativity and enthalpic barriers under isostability conditions. 
!      Phys. Biol. 2 (4), S75-S85 (2005).


 
!  **********************************************************

      PROGRAM MOLECULARDYNAMICS
      
!  *********************************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      character*50 title, date, filename1, filename2, filename
      character*50 term, errmess

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)

      parameter(nbinmax=105)
      parameter(nEbinmax=105)
      parameter(nRbinmax=105)

      dimension xinit(MAXN),yinit(MAXN),zinit(MAXN)
      dimension xsave(MAXN),ysave(MAXN),zsave(MAXN)

!   number of particles for chain1, moving chains, chain1&2, and number of unbonding term (i<j)
!  iFlagMov=0 move chain 1 without influece from chain 2;  =1: move chain 1 only;  =2: move both chains
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
	  
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)
      common/velocity/vx(MAXN),vy(MAXN),vz(MAXN)

      common/forcetotalnew/fx(MAXN),fy(MAXN),fz(MAXN)
      common/forcetotalold/fxo(MAXN),fyo(MAXN),fzo(MAXN)

      common/forcestretching/fxr(MAXN),fyr(MAXN),fzr(MAXN)
      common/forcebending/fxth(MAXN),fyth(MAXN),fzth(MAXN)
      common/forcetorsion/fxph(MAXN),fyph(MAXN),fzph(MAXN)
      common/forceunbonded/fxun(MAXN),fyun(MAXN),fzun(MAXN)

      common/forcerandnew/frandx(MAXN),frandy(MAXN),frandz(MAXN)
      common/forcerandold/frandxo(MAXN),frandyo(MAXN),frandzo(MAXN)
      
      common/nativeinfo/rbond_nat(MAXN),runbond_nat(MAXUNBO),theta_nat(MAXN),dihedral_nat(MAXN)

      common/gyrationradius/gyr	  

      common/constants/pi,boltz,avsn0,gamma,amass,sigma_ij
      common/variables/temp,dt,nstep,nadim,nsnap,gm,enscale

      common/contactmap/iun(MAXUNBO),jun(MAXUNBO),kunbond(MAXUNBO),nQnative_f1,nQnative_f2,nQnative_b

      common/randomterms/c_0,c_1,c_2,randconst
      common/interacparam/ck_r,ck_tht,ck_phi1,ck_phi3,epsil1,epsil2,epsil

      common/LagrangeData/ga1_f1, ga2_f1, ga1_f2, ga2_f2, ga1_b, ga2_b, gQ0_f1, gQ0_f2, gQ0_b, gr0, gdr, gQ_f, gQ_f1, gQ_f2, gQ_b, eGr
      common/LagrangeData2/ga1_w, ga2_w, gQ0_w, gQ_w, alpha_Qw

      common/OutputConform/nConformOuput, outQf1_i, outQf1_f, outQf2_i, outQf2_f, iQb_i, iQb_f, nOutput0, ndOutput
!  output nConfrom conformations with outQf1_i<=Qf1<=outQf1_f and outQf2_i<=Qf2<=outQf2_f and iQb_i<=Qb<=iQb_b
!  noutput0, ndOutput: the initial step and gap step of the output


      common/TransCoefSet/nConform,nCon0, nRunConf, gQbnativ, gQbdenatural, gQf1nativ, gQf1denatural, gQf2nativ, gQf2denatural


      common/Setdtgm/ s_dt, s_gm, iFixKr
      common/Solvation/k_sol, n_sol, m_sol, epsilon_p, epsilon_pp, pseudoQ_f,pseudoQ_b, dr_sol, ddr_sol

      common/WithWithoutSolvation/ Is_Solvation
      dimension ioctsd(4)

!  histogram of (Q_f,Q_b), (Q_f), (Q_b)
      common/HistogramData_fb/nbin_f,nbin_b,dbin_f,dbin_b,vbin0,nbinsnap,nbinsnap0
      common/HistogramData_FB/PFBbin(nbinmax,nbinmax), PFbin(nbinmax), PBbin(nbinmax)

!  histogram of (Q_f,Q_b,E)
      common/EQ_fbHistogramData/nEbin,dEbin,vEbin0, IsEbin
      common/EQ_FBHistogramData/PFBEbin(nbinmax,nbinmax,nEbinmax)

!  histogram of (Q_b,E_b)
      common/EbQbHistogramParameter/nEbbin,dEbbin,vEbbin0, IsEbbin
      common/EbQbHistogramData/PQbEbbin(nbinmax,nEbinmax)

!  histogram of (Q_f,Q_b,R)
      common/RQ_fbHistogramData/nRbin,dRbin,vRbin0, IsRbin
      common/RQ_FBHistogramData/PFBRbin(nbinmax,nbinmax,nRbinmax)

!  histogram of (Q_f1, Q_f2, Q_b)
      common/HistogramData_FFB/PFFBbin(nbinmax,nbinmax,nbinmax)

!  histogram of (Q_w), (Q_w, Q_b) for construction of F(Q_b), (Q_w, R, Qb=0?) for F(R) and K, (Q_w, E_b, Qb=0?) for K changing with epsilon_b
      common/Q_wHistogramParam/nWbin, dWbin, vWbin0, IsWbin, cri_Qb
      common/Q_wHistogramData/PQwbin(nbinmax), PQwQbbin(nbinmax,nbinmax), PQwQfbin(nbinmax,nbinmax), &
                  PQwRbin(nbinmax,nRbinmax,2), PQwEbbin(nbinmax,nEbinmax,2)

!  for periodic boundary condiction,R is the distance to the center of the box	  
	  common/pbc/pL,dl,R

!  Rescaling the interaction strength of interchain interaction and intrachain interaction 
      common/bias/Alpha1,Alpha2,Beta,Delta
      common/non_native/CritR_non, gQ_non_f, gQ_non_b



      write(*, '('' Enter input file name: '')')
      read (*, *) filename
      open ( 100, file = filename, status = 'old' )
   
        read ( 100, '( 2a14 )' ) title, date
!        write( *,   '( 2a14 )' ) title, date
!        write( *, '()' )
        read ( 100, *, end = 10 ) filename1
        
!-------------read chains -------------------      
        read ( 100, * ) npart1, nparttol, nResGap, iFlagMov
        if(iFlagMov.eq.1 .or. iFlagMov.eq.0) then
          npartM=npart1
        else
          npartM=nparttol
        endif
        if(npart1.gt.nparttol) call myErr('Err: npart1.gt.nparttol')
        if(npart1.lt.1 .or. npart1.gt.MAXN) call myErr('Err: npart1.lt.1 .or. .gt.MAXN')
        if(nparttol.gt.MAXN) call myErr('Err: nparttol.gt.MAXN')

!  file of the initial conformation: (InitialConf.dat)
        read ( 100, * ) filename2
	    open ( 20, file = filename2, status = 'old' )
!  Native conformation: (NativeConf.dat)
        read ( 100, * ) filename2
        open ( 2, file = filename2, status = 'old' )
!  contact map file: (appNCS.dat)
        read ( 100, * ) filename2
        open ( 4, file = filename2, status = 'old' )
        read ( 100, * ) term

!-------------read parameters-------------------
        read (100,'(4i4)') ioctsd
        
        read (100,*) epsil, epsil1, epsil2, enscale
        read (100,*) ck_r, ck_tht, ck_phi1, ck_phi3
        read (100,*) sigma_ij, amass, gamma
        read (100,*) avsn0, boltz
	    read (100,*) gr0, gdr
        read (100,*) s_dt, s_gm, iFixKr        
        read (100,*) k_sol, n_sol, m_sol 
        read (100,*) epsilon_p, epsilon_pp
        read (100,*) temp
        read ( 100, * ) ga1_f1, ga2_f1, gQ0_f1
        read ( 100, * ) ga1_f2, ga2_f2, gQ0_f2
	    read ( 100, * ) ga1_b, ga2_b, gQ0_b
	    read ( 100, * ) ga1_w, ga2_w, gQ0_w, alpha_Qw
        read (100,*) Is_Solvation
        read (100,*) nConform, nCon0, nRunConf
        read (100,*) gQbnativ, gQbdenatural, gQf1nativ, gQf1denatural, gQf2nativ, gQf2denatural
        read (100,*) nsnap, nstep
!        read (100,*) ndOutput
        read (100, *) nConformOutput, nOutput0, ndOutput
        read (100, *) outQf1_i, outQf1_f, outQf2_i, outQf2_f, outQb_i, outQb_f
        read (100,*) nbinsnap, nbinsnap0
        read (100,*) nbin_f, nbin_b
        read (100,*) dbin_f, dbin_b, vbin0
        
        read (100,*) IsEbin, nEbin, dEbin, vEbin0
        read (100,*) IsEbbin, nEbbin, dEbbin, vEbbin0
        read (100,*) IsRbin, nRbin, dRbin, vRbin0
        read (100,*) IsWbin, nWbin, dWbin, vWbin0, cri_Qb
        
        read (100,*) pL, dl
        read (100,*) Alpha1, Alpha2, Beta
        read (100,*) Delta, CritR_non	  

        if(nbin_f.gt.nbinmax) call myErr('Err: nbin_f.gt.nbinmax')
        if(nbin_b.gt.nbinmax) call myErr('Err: nbin_b.gt.nbinmax')
        if(nEbin.gt.nEbinmax) call myErr('Err: nEbin.gt.nEbinmax')
        if(nRbin.gt.nRbinmax) call myErr('Err: nRbin.gt.nRbinmax')

!-------------open files to save-------------------
        LF1 = len_trim(filename1) 
 !  file to save output information: (Output*)
        LF2 = LF1 + 6
        filename2(1:6  ) = 'Output'
        filename2(7:LF2) = filename1(1:LF1)
        open ( 9, file = filename2(1:LF2), status = 'replace' )
 !  file to save Histogram P(Q_b): (QbHist*)
        LF2 = LF1 + 6
        filename2(1:6  ) = 'QbHist'
        filename2(7:LF2) = filename1(1:LF1)
        if(iFlagMov.gt.0) open (14,file=filename2(1:LF2), status = 'replace' )
 !  file to save Histogram P(Q_f): (QfHist*)
        LF2 = LF1 + 6
        filename2(1:6  ) = 'QfHist'
        filename2(7:LF2) = filename1(1:LF1)
        open (15,file=filename2(1:LF2), status = 'replace' )
 !  file to save Histogram P(Q_f,Q_b): (QfbHist*)
        LF2 = LF1 + 7
        filename2(1:7  ) = 'QfbHist'
        filename2(8:LF2) = filename1(1:LF1)
        if(iFlagMov.gt.0) open (19,file=filename2(1:LF2), status = 'replace' )
!  file to save the (Q,E)-histogram: (QEHist*.21)
        LF2 = LF1 + 6
        filename2(1:6  ) = 'QEHist'
        filename2(7:LF2) = filename1(1:LF1)
        if(IsEbin.eq.1)  open (23,file=filename2(1:LF2), status = 'replace' )
!  file to save the (Q,R)-histogram: (QRHist*.21)
        LF2 = LF1 + 6
        filename2(1:6  ) = 'QRHist'
        filename2(7:LF2) = filename1(1:LF1)
        if(IsRbin.eq.1)  open (24,file=filename2(1:LF2), status = 'replace' )
!  file to save the (Qf1,Qf2,Qb)-histogram: (QffbHist*.21)
        LF2 = LF1 + 8
        filename2(1:8  ) = 'QffbHist'
        filename2(9:LF2) = filename1(1:LF1)
        if(iFlagMov.eq.2)  open (25,file=filename2(1:LF2), status = 'replace' )
!  file to save the output conformation, d*.mol
        LF2 = LF1 + 5
        filename2(1:1) = 'd'
        filename2(2  :LF1+1) = filename1(1:LF1)
        filename2(LF1+2:LF2) = '.mol'
        open (22,file=filename2(1:LF2), status = 'replace' )
!  file to save the information for output conformations, d*.mol.dat
        LF2 = LF1 + 9
        filename2(1:1) = 'd'
        filename2(2  :LF1+1) = filename1(1:LF1)
        filename2(LF1+2:LF2) = '.mol.dat'
        open(27,file=filename2(1:LF2), status = 'replace' )
!  file to save the output conformation coordinate for further simulation, d*.X 
        LF2 = LF1 + 3
        filename2(2  :LF1+1) = filename1(1:LF1)
        filename2(1:1) = 'd'
        filename2(LF1+2:LF2) = '.X'
        open (26,file=filename2(1:LF2), status = 'replace' )
!  file to save the binding/unbinding time and transmission coefficient, d*.tc
        LF2 = LF1 + 4
        filename2(2  :LF1+1) = filename1(1:LF1)
        filename2(1:1) = 'd'
        filename2(LF1+2:LF2) = '.tc'
        open (28,file=filename2(1:LF2), status = 'replace' )
 !  file to save Histogram P(Q_b,E_b): (QbEbHist*)
        LF2 = LF1 + 8
        filename2(1:8  ) = 'QbEbHist'
        filename2(9:LF2) = filename1(1:LF1)
        if(iFlagMov.gt.0 .and. IsEbbin.eq.1) open (29,file=filename2(1:LF2), status = 'replace' )
 !  file to save Histogram P(Q_w): (QwHist*)
        LF2 = LF1 + 6
        filename2(1:6  ) = 'QwHist'
        filename2(7:LF2) = filename1(1:LF1)
        if(iFlagMov.gt.0 .and. IsWbin.eq.1) open (30,file=filename2(1:LF2), status = 'replace' )
 !  file to save Histogram P(Q_w,Q_b): (QwQbHist*)
        LF2 = LF1 + 8
        filename2(1:8  ) = 'QwQbHist'
        filename2(9:LF2) = filename1(1:LF1)
        if(iFlagMov.gt.0 .and. IsWbin.eq.1) open (31,file=filename2(1:LF2), status = 'replace' )
 !  file to save Histogram P(Q_w,R): (QwRHist*)
        LF2 = LF1 + 7
        filename2(1:7  ) = 'QwRHist'
        filename2(8:LF2) = filename1(1:LF1)
        if(iFlagMov.gt.0 .and. IsWbin.eq.1) open (32,file=filename2(1:LF2), status = 'replace' )
 !  file to save Histogram P(Q_w,E_b): (QwEbHist*)
        LF2 = LF1 + 8
        filename2(1:8  ) = 'QwEbHist'
        filename2(9:LF2) = filename1(1:LF1)
        if(iFlagMov.gt.0 .and. IsWbin.eq.1) open (33,file=filename2(1:LF2), status = 'replace' )
 !  file to save Histogram P(Q_w,Q_f1): (QwEbHist*)
        LF2 = LF1 + 8
        filename2(1:8  ) = 'QwQfHist'
        filename2(9:LF2) = filename1(1:LF1)
        if(iFlagMov.gt.0 .and. IsWbin.eq.1) open (34,file=filename2(1:LF2), status = 'replace' )



        write (*,'('' Length of chain 1 (npart1): '',t50, i10)') npart1
        write (*,'('' Number of moving chains (iFlagMov): '',t50, i10)') iFlagMov
        write (*,'('' Total length (nparttol): '',t50, i10)') nparttol
        write (*,'('' Chain-Residue number gap (nResGap): '',t50, i10)') nResGap
        write (*,'('' interaction strength epsilon'',t50, f10.5)') epsil
        write (*,'('' interaction strength epsilon_1'',t50, f10.5)') epsil1
        write (*,'('' interaction strength epsilon_2'',t50, f10.5)') epsil2
        write (*,'('' interaction strength ck_r'',t50, f10.5)') ck_r
        write (*,'('' interaction strength ck_tht'',t50, f10.5)') ck_tht
        write (*,'('' interaction strength ck_phi1'',t50, f10.5)') ck_phi1
        write (*,'('' interaction strength ck_phi3'',t50, f10.5)') ck_phi3
        write (*,'('' interaction strength scaling param. '',t50, f10.5)')enscale
        write (*,'('' Solvation k: '',t50, i10)') k_sol
        write (*,'('' Solvation n: '',t50, i10)') n_sol
        write (*,'('' Solvation m: '',t50, i10)') m_sol
        write (*,'('' Solvation epsilon pri: '',t50, f10.5)') epsilon_p
        write (*,'('' Solvation epsilon pripri:'',t50, f10.5)')epsilon_pp
        write (*,'('' dt:'',t50,f10.6)')  s_dt
        write (*,'('' gamma:'',t50,f10.6)')  s_gm
        write (*,'('' Is Fix Kr:'',t50,i10)')  iFixKr
        
        write (*,'('' temperature'',t50, f10.5)') temp
        write (*,'('' cut-off'',t50, f10.5)') gamma
        write (*,'('' hard core distance'',t50, f10.5)') sigma_ij
        write (*,'('' total simulation step'',t50, i13)') nstep
        write (*,'('' frequency of snapshots '',t50, i10)') nsnap
        write (*,'('' Avogadro number'',t50, f10.5)') avsn0
        write (*,'('' Boltzmann constant'',t50, f10.5)') boltz
        write (*,'('' random seed'',t50, 4i4)') ioctsd
        
        write (*,'('' Lagrange constraint a1_f1: '',t50, f10.5)') ga1_f1
        write (*,'('' Lagrange constraint a2_f1: '',t50, f10.5)') ga2_f1
        write (*,'('' Lagrange constraint Q0_f1: '',t50, f10.5)') gQ0_f1

        write (*,'('' Lagrange constraint a1_f2: '',t50, f10.5)') ga1_f2
        write (*,'('' Lagrange constraint a2_f2: '',t50, f10.5)') ga2_f2
        write (*,'('' Lagrange constraint Q0_f2: '',t50, f10.5)') gQ0_f2
        
	    write (*,'('' Lagrange constraint a1_b: '',t50, f10.5)') ga1_b
        write (*,'('' Lagrange constraint a2_b: '',t50, f10.5)') ga2_b
	    write (*,'('' Lagrange constraint Q0_b: '',t50, f10.5)') gQ0_b

	    write (*,'('' Lagrange constraint a1_b: '',t50, f10.5)') ga1_2
        write (*,'('' Lagrange constraint a2_b: '',t50, f10.5)') ga2_2
	    write (*,'('' Lagrange constraint Q0_b: '',t50, f10.5)') gQ0_2
	    write (*,'('' Parameter for Q_w as alpha_Qw: '',t50, f10.5)') alpha_Qw
        
        write (*,'('' Lagrange constraint r0: '',t50, f10.5)') gr0
        write (*,'('' Lagrange constraint dr: '',t50, f10.5)') gdr
        
        write (*,'('' nConformOutput:'',t50,i10)') nConformOutput
        write (*,'('' nOutput0 :'',t50,i10)') nOutput0
        write (*,'('' ndOutput:'',t50,i10)') ndOutput
        write (*,'('' outQf1_i: '',t50, f10.5)') outQf1_i
        write (*,'('' outQf1_f: '',t50, f10.5)') outQf1_f
        write (*,'('' outQf2_i: '',t50, f10.5)') outQf2_i
        write (*,'('' outQf2_f: '',t50, f10.5)') outQf2_f
        write (*,'('' outQb_i: '',t50, f10.5)') outQb_i
        write (*,'('' outQb_f: '',t50, f10.5)') outQb_f
        
        write (*,'('' Histogram nbin_f: '',t50, i10)') nbin_f
        write (*,'('' Histogram dbin_f: '',t50, f10.5)') dbin_f
	    write (*,'('' Histogram nbin_b: '',t50, i10)') nbin_b
        write (*,'('' Histogram dbin_b: '',t50, f10.5)') dbin_b
        
        write (*,'('' num. of conformation:'',t50,i10)')nConform 
        write (*,'('' num. of conformation to skip:'',t50,i10)')nCon0 
        write (*,'('' nRunConf'',t50,i10)') nRunConf
        write (*,'('' gQbnativ'',t50,f10.3)') gQbnativ
        write (*,'('' gQbdenatural:'',t50,f10.3)')  gQbdenatural
        write (*,'('' gQf1nativ'',t50,f10.3)') gQf1nativ
        write (*,'('' gQf1denatural:'',t50,f10.3)')  gQf1denatural
        write (*,'('' gQf2nativ'',t50,f10.3)') gQf2nativ
        write (*,'('' gQf2denatural:'',t50,f10.3)')  gQf2denatural


        write (*,'('' Force scheme: '',t50,i10)') Is_Solvation
        write (*,'('' If save E-Histogram: '',t50,i10)') IsEbin
        write (*,'('' E-Histogram nbin: '',t50, i10)') nEbin
        write (*,'('' E-Histogram dbin: '',t50, f10.5)') dEbin
        write (*,'('' E-Histogram vbin0: '',t50, f10.5)') vEbin0
        write (*,'('' If save R-Histogram: '',t50,i10)') IsRbin
        write (*,'('' R-Histogram nbin: '',t50, i10)') nRbin
        write (*,'('' R-Histogram dbin: '',t50, f10.5)') dRbin
        write (*,'('' R-Histogram vbin0: '',t50, f10.5)') vRbin0
        write (*,'('' If save Qw-Histogram: '',t50,i10)') IsWbin
        write (*,'('' Qw-Histogram nbin: '',t50, i10)') nWbin
        write (*,'('' Qw-Histogram dbin: '',t50, f10.5)') dWbin
        write (*,'('' Qw-Histogram vbin0: '',t50, f10.5)') vWbin0
        write (*,'('' cri_Qb for Qw-Histogram: '',t50, f10.5)') cri_Qb
	    write (*,'('' periodic boundary condiction size:'',t50, f10.5)') pL
	    write (*,'('' periodic boundary condiction buffer:'',t50, f10.5)') dl
        write (*,'('' Ratio of intrachain interaction 1:'',t50, f10.5)') Alpha1
        write (*,'('' Ratio of intrachain interaction 2:'',t50, f10.5)') Alpha2
	    write (*,'('' Ratio of interchain interaction:'',t50, f10.5)') Beta
	    write (*,'('' Delta:'',t50, f10.5)') Delta
	    write (*,'('' CritR_non:'',t50, f10.5)') CritR_non



        read (2,*)(x(j),y(j),z(j),j=1,nparttol)

        nunbond=0
        do while(nunbond.lt.MAXUNBO)
          read(4,*, ERR=17, END=17) item1, item2, item3
		  if(item1.gt.npart1) then
		    if(item1.lt.npart1+nResGap) call myErr('contact map error 4.')
		    item1=item1-nResGap
		  endif
		  if(item2.gt.npart1) then
		    if(item2.lt.npart1+nResGap) call myErr('contact map error 5.')
		    item2=item2-nResGap
		  endif		   
          if(item1.le.0.or.item1.gt.nparttol) call myErr('contact map error 1.')
          if(item2.le.0.or.item2.gt.nparttol) call myErr('contact map error 2.')
          if(item1.eq.item2) call myErr('contact map error 3.')
          nunbond=nunbond+1
          if(item1.lt.item2) then
            iun(nunbond)=item1
            jun(nunbond)=item2
          else
            iun(nunbond)=item2
            jun(nunbond)=item1
          endif
          kunbond(nunbond)=item3
          if(item1.gt.npartM .and. item2.gt.npartM .or. &
              iFlagMov.eq.0 .and. (item1.gt.npart1 .or. item2.gt.npart1)) nunbond=nunbond-1
        end do
 17     continue
        write (*,'('' Effective Number of unbond: '',i7)') nunbond
!	    read (4,*)( iun(j), jun(j), kunbond(j), runbond_nat(j), j = 1, nunbond )
        
        call setrn(ioctsd)
        call nativeinformation
        call intpar(enerkin)

        write (*,'('' Native Qf:'',i4,i4,'';   Qb:'',i4)') nQnative_f1,nQnative_f2,nQnative_b


!  total number of folding/unfolding events
        nFoldTol=0
        nunFoldTol=0
        aveNadim=0
        
        do 316 j=1,nCon0
          do 315 i=1,nparttol
            read(20,*) xinit(i),yinit(i),zinit(i)
 315      continue
 316    continue


!===================== simulation start =====================
        do 2000 nConfcount=nCon0+1, nConform
          do 320 i=1,nparttol
            read(20,*) xinit(i),yinit(i),zinit(i)
 320      continue

          nFoldSub=0
          nunFoldSub=0
        
          do 1500 nRuncount=1, nRunConf
            do 330 i=1, nparttol
              x(i)=xinit(i)
              y(i)=yinit(i)
              z(i)=zinit(i)
 330        continue
           
            call origin_adjust
            call InitVel(enerkin)
            call force(e_pot,e_unbond_tot,e_bind_tot,e_tors_tot,e_bend_tot,e_bond_tot)
            call RANTERM      
        
            do 319 kl=1,npartM
              fxo(kl)=fx(kl) 
              fyo(kl)=fy(kl) 
              fzo(kl)=fz(kl) 
              frandxo(kl)=frandx(kl)
              frandyo(kl)=frandy(kl)
              frandzo(kl)=frandz(kl)
  319       continue                                                   
        
            nOutputCount=0
            nOutputGap=-1
            nadim=-1
          

! main dynamics cycle
            do 100 nadim_new=0,nstep
              nadim=nadim+1
              call verlet(enerkin,e_pot)
            
            
!-------------save data for purpose of restore-------------
              if(mod(nadim,nsnap).eq.0 .and. vx(1).ge.-1e6 .and. vx(1).le.1e6) then
                do 40 kl=1,nparttol
                  xsave(kl)=x(kl) 
                  ysave(kl)=y(kl) 
                  zsave(kl)=z(kl) 
 40             continue
                  
                nadim_old=nadim
                nOutGap_old=nOutputGap
              endif
!  restore
              if(.not.(vx(1).ge.-1e6 .and. vx(1).le.1e6)) then
                do 45 i=1, nparttol
                  x(i)=xsave(i)
                  y(i)=ysave(i)
                  z(i)=zsave(i)
 45             continue
             
		        call InitVel(enerkin)
                call force(e_pot,e_unbond_tot,e_bind_tot,e_tors_tot,e_bend_tot,e_bond_tot)
                call RANTERM      
                do 48 kl=1,npartM
                  fxo(kl)=fx(kl) 
                  fyo(kl)=fy(kl) 
                  fzo(kl)=fz(kl) 
                  frandxo(kl)=frandx(kl)
                  frandyo(kl)=frandy(kl)
                  frandzo(kl)=frandz(kl)
 48             continue                                                   
                nadim=nadim_old
                nOutputGap=nOutGap_old
              endif
!-------------save data for purpose of restore end----------          
                 

!--------------output conformation-------------
              nOutputGap=nOutputGap+1

              if ( nadim == 0 ) write ( 27, '( '' nOutputCount,  nadim,   gQ_f,     gQ_b,     R'' )' )
              if(nOutputCount.lt.nConformOutput .and. nadim.ge.nOutput0 .and. nOutputGap.ge.ndOutput .and. &
                   gQ_f1.ge.outQf1_i .and.  gQ_f1.le.outQf1_f .and. &
                   gQ_f2.ge.outQf2_i .and.  gQ_f2.le.outQf2_f .and. &
                   gQ_b.ge.outQb_i .and.  gQ_b.le.outQb_f) then
                nOutputGap=0
                call write_mol(nOutputCount)
                write ( 27, '( i7, i13, 3f10.3 )' ) nOutputCount, nadim, gQ_f, gQ_b, R
                do 50 i=1, nparttol
                  write(26, *) x(i),y(i),z(i)
 50             continue
                write(26,'('' '')')
                write(*,'(''n='',i3,'', nadim='',i8,'', gQ='',f7.3, '', old gQ='', 1i3)')nOutputCount,nadim,gQ,natcont
              endif
!--------------output conformation end---------

            if(gQ_f1.le.gQf1denatural .and. gQ_f2.le.gQf2denatural .and. gQ_b.le.gQbdenatural & 
			   .or. gQ_f1.ge.gQf1nativ .and. gQ_f2.ge.gQf2nativ .and. gQ_b.ge.gQbnativ) then
              goto 1000
            endif

 100        continue 
 1000     continue
          aveNadim=aveNadim+nadim
          if(.not.(vx(1).ge.-1e8 .and. vx(1).le.1e8)) then
            k=0
          else if(gQ_f1.le.gQf1denatural .and. gQ_f2.le.gQf2denatural .and. gQ_b.le.gQbdenatural) then
            k=-1
            nunFoldSub=nunFoldSub+1
            nunFoldTol=nunFoldTol+1
          else if(gQ_f1.ge.gQf1nativ .and. gQ_f2.ge.gQf2nativ .and. gQ_b.ge.gQbnativ) then
            k=1
            nFoldSub=nFoldSub+1
            nFoldTol=nFoldTol+1
          else
            k=0
          endif
          write(28,'(3i5,i13,3f10.3)') nConfcount, nRuncount, k, nadim, gQ_b, gQ_f1, gQ_f2
 1500     continue
          write(28,'(''# Conform No.:'',i5,'',   Folding times:'',i5,'',   Transmission coefficient:'',f10.3,'',   unFolding times:'',i5)')  &
		            nConfcount,nFoldSub,1.0*nFoldSub/nRunConf, nunFoldSub
 2000   continue
	    aveNadim=aveNadim/((nConform-nCon0)*nRunConf)
        write(28,'(''# Total folding times:'',i5,'',   Total unfolding times:'',i5)') nFoldTol, nunFoldTol
	    write(28,'(''# Average Running Time :'',f14.3)') aveNadim


!===================== simulation end =====================


! output target
        call write_target(nOutputCount)

! output histogram
        call write_histogram

	    close(2)
	    close(4)
	    close(20)
	    close(5)
	    close(9)
	    if(iFlagMov.gt.0) close(14)
	    close(15)
	    if(iFlagMov.gt.0) close(19)
	    close(21)
	    close(22)
	    if(IsEbin.eq.1) close(23)
	    if(IsRbin.eq.1) close(24)
	    if(iFlagMov.eq.2) close(25)
        close(26)
        close(27)
        close(28)
	    if(iFlagMov.gt.0 .and. IsWbin.eq.1) then
          close(30)
          close(31)
          close(32)
          close(33)
          close(34)
	    endif
	  

 
10    continue
      
      
      end


!  *****************************************************
      SUBROUTINE myErr(errmess)  
!  ****************************************************
      implicit real*8 (a-h,o-z) 
      implicit integer*4 (i-n)

      character errmess(100)

      write(*,*) errmess
      write(9,*) errmess
      call exit(0)
      
      return
      end



! *******************************************

      SUBROUTINE write_histogram

! *******************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      
	  parameter(nbinmax=105)
      parameter(nEbinmax=105)
	  parameter(nRbinmax=105)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond

!  histogram of (Q_f,Q_b)
      common/HistogramData_fb/nbin_f,nbin_b,dbin_f,dbin_b,vbin0,nbinsnap,nbinsnap0
      common/HistogramData_FB/PFBbin(nbinmax,nbinmax), PFbin(nbinmax), PBbin(nbinmax)

!  histogram of (Q_f,Q_b,E)
      common/EQ_fbHistogramData/nEbin,dEbin,vEbin0, IsEbin
      common/EQ_FBHistogramData/PFBEbin(nbinmax,nbinmax,nEbinmax)

!  histogram of (Q_b,E_b)
      common/EbQbHistogramParameter/nEbbin,dEbbin,vEbbin0, IsEbbin
      common/EbQbHistogramData/PQbEbbin(nbinmax,nEbinmax)

!  histogram of (Q_f,Q_b,R)
      common/RQ_fbHistogramData/nRbin,dRbin,vRbin0, IsRbin
      common/RQ_FBHistogramData/PFBRbin(nbinmax,nbinmax,nRbinmax)

!  histogram of (Q_f1, Q_f2, Q_b)
      common/HistogramData_FFB/PFFBbin(nbinmax,nbinmax,nbinmax)

!  histogram of (Q_w), (Q_w, Q_b) for construction of F(Q_b), (Q_w, R, Qb=0?) for F(R) and K, (Q_w, E_b, Qb=0?) for K changing with epsilon_b
      common/Q_wHistogramParam/nWbin, dWbin, vWbin0, IsWbin, cri_Qb
      common/Q_wHistogramData/PQwbin(nbinmax), PQwQbbin(nbinmax,nbinmax), PQwQfbin(nbinmax,nbinmax), &
                  PQwRbin(nbinmax,nRbinmax,2), PQwEbbin(nbinmax,nEbinmax,2)


!  histogram result saving to file

	  do 20 i=1,nbin_f
        write(15, '(f6.2, f12.1)') vbin0+(i-0.5)*dbin_f, PFbin(i)
 20   continue
 
 
      if (iFlagMov.le.0) goto 25

	  do 1001 i=1,nbin_b
        sumQ_b=0
	    do 1002 j=1,nbin_f
		  sumQ_b=sumQ_b+PFBbin(j,i)
 1002   continue
        write(14, '(f6.2, f12.1)') vbin0+(i-0.5)*dbin_b,sumQ_b   
 1001 continue
 25   continue

     
	  do 1003 i=1, nbin_f
        do 1004 j=1, nbin_b        
          if(iFlagMov.gt.0 .and. PFBbin(i,j).gt.1e-5)  write(19, '(2f7.2, f12.1)') vbin0+(i-0.5)*dbin_f,vbin0+(j-0.5)*dbin_b,PFBbin(i,j)
          if ( IsEbin == 0 ) goto 1007
		  do 1005 k=1,nEbin
		    if(PFBEbin(i,j,k).gt.1e-5) then
		      write(23, '(2f7.2, f8.2, f12.1)') vbin0+(i-0.5)*dbin_f,vbin0+(j-0.5)*dbin_b,vEbin0+(k-0.5)*dEbin,PFBEbin(i,j,k)
		    endif
 1005     continue
 1007     continue

          if ( IsRbin == 0 ) goto 1008
          do 1006 l=1,nRbin 
		    if(PFBRbin(i,j,l).gt.1e-5) then
		      write(24, '(3f7.2, f12.1)') vbin0+(i-0.5)*dbin_f,vbin0+(j-0.5)*dbin_b,vRbin0+(l-0.5)*dRbin,PFBRbin(i,j,l)
		    endif
 1006     continue
 1008     continue
 
          if ( iFlagMov .ne. 2 ) goto 1010
          do 1009 l=1,nbin_f
		    if(PFFBbin(i,l,j).gt.1e-5) then
		      write(25, '(3f7.2, f12.1)') vbin0+(i-0.5)*dbin_f,vbin0+(l-0.5)*dbin_f,vbin0+(j-0.5)*dbin_b,PFFBbin(i,l,j)
		    endif
 1009     continue
 1010     continue
 
 1004   continue
 1003 continue


          if ( iFlagMov.le.0 .or. IsEbbin.eq.0 ) goto 40
          do 35 j=1, nbin_b        
		  do 30 k=1,nEbbin
		    if(PQbEbbin(j,k).gt.1e-5) then
		      write(29, '(2f8.2, f12.1)') vbin0+(j-0.5)*dbin_b, vEbbin0+(k-0.5)*dEbbin,PQbEbbin(j,k)
		    endif
 30       continue
 35       continue
 40       continue


          if ( iFlagMov.le.0 .or. IsWbin.eq.0 ) goto 180
          do 150 i=1, nWbin        
		    write(30, '(f8.2, f12.1)') vWbin0+(i-0.5)*dWbin, PQwbin(i)
            do 110 j=1, nbin_b        
		      if(PQwQbbin(i,j).gt.1e-5) then
		        write(31, '(2f8.2, f12.1)') vWbin0+(i-0.5)*dWbin, vbin0+(j-0.5)*dbin_b, PQwQbbin(i,j)
		      endif
 110        continue
            do 120 l=1,nRbin 
		      if(PQwRbin(i,l,1).gt.1e-5 .or. PQwRbin(i,l,2).gt.1e-5) then
		        write(32, '(2f8.2, 2f12.1)') vWbin0+(i-0.5)*dWbin, vRbin0+(l-0.5)*dRbin, PQwRbin(i,l,1), PQwRbin(i,l,2)
		      endif
 120        continue
		    do 130 k=1,nEbbin
		      if(PQwEbbin(i,k,1).gt.1e-5 .or. PQwEbbin(i,k,2).gt.1e-5) then
		        write(33, '(2f8.2, 2f12.1)') vWbin0+(i-0.5)*dWbin, vEbbin0+(k-0.5)*dEbbin, PQwEbbin(i,k,1), PQwEbbin(i,k,2)
		      endif
 130        continue
            do 140 j=1, nbin_f        
		      if(PQwQfbin(i,j).gt.1e-5) then
		        write(34, '(2f8.2, f12.1)') vWbin0+(i-0.5)*dWbin, vbin0+(j-0.5)*dbin_f, PQwQfbin(i,j)
		      endif
 140        continue
 150      continue
 180      continue


      end



! *******************************************

      SUBROUTINE write_mol(nOutputCount)

! output conformations in MOL format
! *******************************************

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)


      nOutputCount = nOutputCount + 1
      
	  write(22,'(''Molecule-'',i7)')	 nOutputCount
	  write(22,'(''  ViewerPro         3D                             0'',/)')
	  if(iFlagMov.eq.1 .or. iFlagMov.eq.0) then
	    write(22,'(2i3,''  0  0  0  0  0  0  0  0999 V2000'')') npart1,npart1-1
	  else
	  	write(22,'(2i3,''  0  0  0  0  0  0  0  0999 V2000'')') nparttol,nparttol-2
	  endif
	  
      do 50 i=1, npartM
        write(22, '(3f10.4,'' C   0  0  0  0  0  0  0  0  0  0'')') x(i),y(i),z(i)
 50   continue
	  do 211 i=1,npartM-1
	    j=i+1
        if(j.le.npart1.or.i.gt.npart1) write(22,'(2i3,''  1  0  0  0'')')i,j
211   continue
      write(22,'(''M  END'')')
	  write(22,'(''$$$$'')')


      return
      
      end



! *******************************************

      SUBROUTINE write_target(nOutputCount)

! output conformations in MOL format
! *******************************************

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)


       nOutputCount=nOutputCount+1
	   write(22,'(''Molecule-'',i7)')	 nOutputCount
	   write(22,'(''  ViewerPro         3D                             0'',/)')
	   write(22,'(2i3,''  0  0  0  0  0  0  0  0999 V2000'')') nparttol-npart1,nparttol-npart1-1
       do 511 i=npart1+1, nparttol
         write(22, '(3f10.4,'' C   0  0  0  0  0  0  0  0  0  0'')') x(i),y(i),z(i)
 511   continue
	   do 512 i=1,nparttol-npart1-1
	     j=i+1
         write(22,'(2i3,''  1  0  0  0'')')i,j
512    continue
       write(22,'(''M  END'')')
	   write(22,'(''$$$$'')')

      return
      
      end





! *******************************************

      SUBROUTINE origin_adjust

! place the target to (0,0,0)
! *******************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)

      if(nparttol-npart1.le.0) goto 300

      x00=0.0
      y00=0.0
      z00=0.0

      do 281 j=npart1+1,nparttol
        x00=x00+xparticle(j)
        y00=y00+y(j)
        z00=z00+z(j)
  281 continue

      x00=x00/(nparttol-npart1)
      y00=y00/(nparttol-npart1)
      z00=z00/(nparttol-npart1)

	  do 282 j=1,nparttol
	    x(j)=x(j)-x00
		y(j)=y(j)-y00
		z(j)=z(j)-z00
  282 continue

  300 continue
      
      end



      
!  *******************************************
      SUBROUTINE INTPAR(enerkin)
!  *******************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      
      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
	  parameter(nbinmax=105)
      parameter(nEbinmax=105)
	  parameter(nRbinmax=105)

      common/velocity/vx(MAXN),vy(MAXN),vz(MAXN)
      common/constants/pi,boltz,avsn0,gamma,amass,sigma_ij
      common/variables/temp,dt,nstep,nadim,nsnap,gm,enscale
      common/randomterms/c_0,c_1,c_2,randconst
      common/interacparam/ck_r,ck_tht,ck_phi1,ck_phi3,epsil1,epsil2,epsil
      common/Setdtgm/ s_dt, s_gm, iFixKr
      common/Solvation/k_sol, n_sol, m_sol, epsilon_p, epsilon_pp, pseudoQ_f,pseudoQ_b, dr_sol, ddr_sol

!  histogram of (Q_f,Q_b)
      common/HistogramData_fb/nbin_f,nbin_b,dbin_f,dbin_b,vbin0,nbinsnap,nbinsnap0
      common/HistogramData_FB/PFBbin(nbinmax,nbinmax), PFbin(nbinmax), PBbin(nbinmax)

!  histogram of (Q_f,Q_b,E)
      common/EQ_fbHistogramData/nEbin,dEbin,vEbin0, IsEbin
      common/EQ_FBHistogramData/PFBEbin(nbinmax,nbinmax,nEbinmax)

!  histogram of (Q_b,E_b)
      common/EbQbHistogramParameter/nEbbin,dEbbin,vEbbin0, IsEbbin
      common/EbQbHistogramData/PQbEbbin(nbinmax,nEbinmax)

!  histogram of (Q_f,Q_b,R)
      common/RQ_fbHistogramData/nRbin,dRbin,vRbin0, IsRbin
      common/RQ_FBHistogramData/PFBRbin(nbinmax,nbinmax,nRbinmax)

!  histogram of (Q_f1, Q_f2, Q_b)
      common/HistogramData_FFB/PFFBbin(nbinmax,nbinmax,nbinmax)

!  histogram of (Q_w), (Q_w, Q_b) for construction of F(Q_b), (Q_w, R, Qb=0?) for F(R) and K, (Q_w, E_b, Qb=0?) for K changing with epsilon_b
      common/Q_wHistogramParam/nWbin, dWbin, vWbin0, IsWbin, cri_Qb
      common/Q_wHistogramData/PQwbin(nbinmax), PQwQbbin(nbinmax,nbinmax), PQwQfbin(nbinmax,nbinmax), &
                  PQwRbin(nbinmax,nRbinmax,2), PQwEbbin(nbinmax,nEbinmax,2)

      common/non_native/CritR_non, gQ_non_f, gQ_non_b


      temp=temp*epsil/boltz
      timeunit=sqrt(amass/epsil)*sigma_ij
!       dt=0.005*timeunit
!       gm=0.05/timeunit
      dt=s_dt*timeunit
      gm=s_gm/timeunit

      c_0=dt*gm/(2.0*amass)
      c_1=(1.0-c_0)*(1.0-c_0+c_0**2)
      c_2=(dt/(2*amass))*(1.0-c_0+c_0**2)
      randconst=sqrt(2.0*boltz*temp*gm/dt)
 
!  no change ck_r    -ZRLiu
!       ck_r   = enscale*ck_r
      if(iFixKr.eq.0) ck_r=enscale*ck_r
      ck_tht = enscale*ck_tht
      ck_phi1= enscale*ck_phi1
      ck_phi3= enscale*ck_phi3
      epsil1 = enscale*epsil1
      epsil2 = enscale*epsil2
      epsilon_p=enscale*epsilon_p
      epsilon_pp=enscale*epsilon_pp

! ccccccccccccccccccccccccccccccccc
!     initiliaze velocities
! ccccccccccccccccccccccccccccccccc

      sumvx=0.0 
      sumvy=0.0 
      sumvz=0.0 

      vm=sqrt(temp*boltz/amass)
      
	  do 1001 i=1,npartM
      vx(i)=vm*gauss(xsi)
      vy(i)=vm*gauss(xsi)
      vz(i)=vm*gauss(xsi)
      sumvx=sumvx+vx(i) 
      sumvy=sumvy+vy(i) 
      sumvz=sumvz+vz(i) 
 1001 continue
      
	  sumvx=sumvx/npartM
      sumvy=sumvy/npartM
      sumvz=sumvz/npartM
      enerkin=0.0
      
	  do 3 i=1,npartM
      vx(i)=vx(i)-sumvx 
      vy(i)=vy(i)-sumvy
      vz(i)=vz(i)-sumvz
      enerkin=enerkin+vx(i)**2+vy(i)**2+vz(i)**2
   3  continue
      
	  enerkin =  0.5 * enerkin * amass 
      tempins=2.0*enerkin/(3.0*npartM*boltz) 

      
!  initialize histogram data
      do 30 i=1,nbin_f
	    PFbin(i)=0.0
	    do 40 j=1,nbin_b
          PBbin(j)=0.0
          PFBbin(i,j)=0.0
		  do 50 k=1,nEbin
		    PFBEbin(i,j,k)=0.0
  50      continue
		  do 55 k=1,nEbbin
		    PQbEbbin(j,k)=0.0
  55      continue
          do 60 l=1,nRbin 
		    PFBRbin(i,j,l)=0.0
  60      continue		   
          do 70 k=1,nbin_f 
		    PFFBbin(i,k,j)=0.0
  70      continue		   
  40    continue
  30  continue  

          if ( iFlagMov.le.0 .or. IsWbin.eq.0 ) goto 180
          do 150 i=1, nWbin        
		    PQwbin(i)=0
            do 110 j=1, nbin_b        
		        PQwQbbin(i,j)=0
 110        continue
            do 120 l=1,nRbin 
              PQwRbin(i,l,1)=0
              PQwRbin(i,l,2)=0
 120        continue
		    do 130 k=1,nEbbin
		        PQwEbbin(i,k,1)=0
		        PQwEbbin(i,k,2)=0
 130        continue
            do 140 j=1, nbin_f       
		        PQwQfbin(i,j)=0
 140        continue
 150      continue
 180      continue

      return
      end


!  ******************************************
      SUBROUTINE FORCE(e_pot,e_unbond_tot,e_bind_tot,e_tors_tot,e_bend_tot,e_bond_tot)
!  *****************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
 
      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/forcetotalnew/fx(MAXN),fy(MAXN),fz(MAXN)
      common/forcestretching/fxr(MAXN),fyr(MAXN),fzr(MAXN)
      common/forcebending/fxth(MAXN),fyth(MAXN),fzth(MAXN)
      common/forcetorsion/fxph(MAXN),fyph(MAXN),fzph(MAXN)
      common/forceunbonded/fxun(MAXN),fyun(MAXN),fzun(MAXN)

      e_pot=0.0

      do 2  i=1,npartM
      fx(i)=0
      fy(i)=0
      az(i)=0
      fxr(i)=0
      fyr(i)=0
      fzr(i)=0
      fxth(i)=0
      fyth(i)=0
      fzth(i)=0
      fxph(i)=0
      fyph(i)=0
      fzph(i)=0
      fxun(i)=0
      fyun(i)=0
      fzun(i)=0
   2  continue

      e_bond_tot=0.0
      e_bend_tot=0.0
      e_tors_tot=0.0
      e_unbond_tot=0.0
      e_bind_tot=0.0

      call fbond(e_bond_tot)
      call fbend(e_bend_tot)
      call ftorsion(e_tors_tot)
      call funbond(e_unbond_tot,e_bind_tot)

      do 60 i=1,npartM
      fx(i)=fx(i)+fxr(i)+fxth(i)+fxph(i)+fxun(i)
      fy(i)=fy(i)+fyr(i)+fyth(i)+fyph(i)+fyun(i)
      fz(i)=fz(i)+fzr(i)+fzth(i)+fzph(i)+fzun(i)
  60  continue
      e_pot=e_bond_tot+e_bend_tot+e_tors_tot+e_unbond_tot
      return 
      end



!  *****************************************************
      SUBROUTINE FUNBOND(e_unbond_tot,e_bind_tot)
!  *****************************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      common/WithWithoutSolvation/ Is_Solvation

      if(Is_Solvation.eq.0) then
        call FUNBOND_without(e_unbond_tot,e_bind_tot)
      else if(Is_Solvation.eq.1) then
        call FUNBOND_with(e_unbond_tot,e_bind_tot)
      endif

      return
      end




!  *****************************************************
      SUBROUTINE FUNBOND_without(e_unbond_tot,e_bind_tot)
!  *****************************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)
      common/forceunbonded/fxun(MAXN),fyun(MAXN),fzun(MAXN)
      common/nativeinfo/rbond_nat(MAXN),runbond_nat(MAXUNBO),theta_nat(MAXN),dihedral_nat(MAXN)
      common/constants/pi,boltz,avsn0,gamma,amass,sigma_ij
      common/variables/temp,dt,nstep,nadim,nsnap,gm,enscale
      common/contactmap/iun(MAXUNBO),jun(MAXUNBO),kunbond(MAXUNBO),nQnative_f1,nQnative_f2,nQnative_b
      common/interacparam/ck_r,ck_tht,ck_phi1,ck_phi3,epsil1,epsil2,epsil
      common/LagrangeData/ga1_f1, ga2_f1, ga1_f2, ga2_f2, ga1_b, ga2_b, gQ0_f1, gQ0_f2, gQ0_b, gr0, gdr, gQ_f, gQ_f1, gQ_f2, gQ_b, eGr
      common/LagrangeData2/ga1_w, ga2_w, gQ0_w, gQ_w, alpha_Qw

      common/bias/Alpha1,Alpha2,Beta,Delta
      common/non_native/CritR_non, gQ_non_f, gQ_non_b
      
      
      gQ_f=0.0
      gQ_f1=0.0
      gQ_f2=0.0
      gQ_b=0.0
      gQ_w=0.0
      sum_rij=0.0
      gQ_non_f=0.0
      gQ_non_b=0.0
      
      do 40 k=1,nunbond
        i=iun(k)
        j=jun(k)

        e_unbond=0
        e_unbond_drv=0

        xij=(x(i)-x(j))      
        yij=(y(i)-y(j))      
        zij=(z(i)-z(j))
        rij=sqrt(xij**2+yij**2+zij**2)

! cancel for Q_w
!      if ( kunbond(k) == 1 .and. rij >= (2*runbond_nat(k)+3.0) .or. &

       if ( kunbond(k) == 2 .and. rij >= (2*CritR_non+3.0     ) .or. & 
           kunbond(k) == 0 .and. rij >= (1.8*sigma_ij        ) .or. &
           i>npartM .and. j>nparM )    goto 60

      if(kunbond(k) == 1) then
         if((rij/runbond_nat(k)-gr0) .lt. 8*gdr) then
            tem=1/(1+exp((rij/runbond_nat(k)-gr0)/gdr))
            if(i.le.npart1 .and. j.le.npart1) then
              gQ_f1=gQ_f1+tem
            else if(i.le.npart1 .and. j.gt.npart1) then
              gQ_b=gQ_b+tem
            else if(iFlagMov.eq.2) then
              gQ_f2=gQ_f2+tem
            endif
         endif
         if(i.le.npart1 .and. j.gt.npart1) sum_rij=sum_rij+(rij-runbond_nat(k))
         rij10=(runbond_nat(k)/rij)**10.0
         rij12=(runbond_nat(k)/rij)**12.0
         if(i.le.npart1 .and. j.le.npart1) then
           e_unbond = Alpha1 * epsil1 * ( 5 * rij12 - 6 * rij10 )
           e_unbond_drv = Alpha1 * 60 * epsil1 * ( rij12 - rij10 ) / rij
         else if(i.le.npart1 .and. j.gt.npart1) then
           e_unbond = Beta * epsil1 * ( 5 * rij12 - 6 * rij10 )
           e_unbond_drv = Beta * 60 * epsil1 * ( rij12 - rij10 ) / rij
         else if(iFlagMov.eq.2) then
           e_unbond = Alpha2 * epsil1 * ( 5 * rij12 - 6 * rij10 )
           e_unbond_drv = Alpha2 * 60 * epsil1 * ( rij12 - rij10 ) / rij
          endif     
      else if(kunbond(k) == 2) then
         if((rij/CritR_non-gr0) .lt. 8*gdr) then
            tem=1/(1+exp((rij/CritR_non-gr0)/gdr))
            if(i.le.npart1 .and. j.gt.npart1) then
              gQ_non_b=gQ_non_b+tem
            else if(iFlagMov.eq.2 .or. j.le.npart1) then
              gQ_non_f=gQ_non_f+tem
            endif
         endif      	
      	
         rij12=(sigma_ij/rij)**12.0
         rHP  = exp( -1 * (rij - CritR_non )**2 / 2 )  
         if(i.le.npart1 .and. j.le.npart1) then
           Atem = Alpha1
         else if(i.le.npart1 .and. j.gt.npart1) then
!           Atem = sqrt(Alpha1*Alpha2)
           Atem = Beta
         else
           Atem = Alpha2
         endif     
         e_unbond = Atem * epsil2 * rij12 - Delta * epsil1 * rHP
         e_unbond_drv = Atem * 12 * epsil2 * rij12 / rij - Delta * epsil1 * rHP * ( rij - CritR_non )    

!         print *, e_unbond, e_unbond_drv 
      
      else
         rij12=(sigma_ij/rij)**12.0
         if(i.le.npart1 .and. j.le.npart1) then
           Atem = Alpha1
         else if(i.le.npart1 .and. j.gt.npart1) then
!           Atem = sqrt(Alpha1*Alpha2)
           Atem = Beta
         else
           Atem = Alpha2
         endif     
         e_unbond = Atem * epsil2 * rij12
         e_unbond_drv = Atem * 12 * epsil2 * rij12 / rij
      endif
      
      e_unbond_tot=e_unbond_tot+e_unbond
      if(i.le.npart1 .and. j.gt.npart1) e_bind_tot=e_bind_tot+e_unbond

      if(i.le.npartM) then
        fxun(i)=fxun(i)+e_unbond_drv*xij/rij
        fyun(i)=fyun(i)+e_unbond_drv*yij/rij
        fzun(i)=fzun(i)+e_unbond_drv*zij/rij
      endif
      if(j.le.npartM) then
        fxun(j)=fxun(j)-e_unbond_drv*xij/rij
        fyun(j)=fyun(j)-e_unbond_drv*yij/rij
        fzun(j)=fzun(j)-e_unbond_drv*zij/rij
      endif


!      Non-native interaction used:
!      Theoretical and experimental demonstration of
!      the importance of specific nonnative interactions
!      in protein folding
!      PNAS July 22, 2008 vol. 105: 9999¨C10004

   60 continue
   40 continue
   
      if(iFlagMov.ne.2) gQ_f2=0
      gQ_f=gQ_f1+gQ_f2



!  Lagrange constrant potential

      gQ_w=gQ_b-alpha_Qw*sum_rij

      eGr_f1 = ga1_f1 * ( gQ_f1 - gQ0_f1 ) + ga2_f1 * ( gQ_f1 - gQ0_f1 ) ** 2
      if(iFlagMov.eq.2) then
        eGr_f2 = ga1_f2 * ( gQ_f2 - gQ0_f2 ) + ga2_f2 * ( gQ_f2 - gQ0_f2 ) ** 2
      else
        eGr_f2 = 0.0
      endif
      eGr_b = ga1_b * ( gQ_b - gQ0_b ) + ga2_b * ( gQ_b - gQ0_b ) ** 2
      eGr_f=eGr_f1+eGr_f2
      eGr_w = ga1_w * ( gQ_w - gQ0_w ) + ga2_w * ( gQ_w - gQ0_w ) ** 2
      eGr   = eGr_f + eGr_b +eGr_w

      if(eGr.lt.1e-5 .and. eGr.gt.-1e-5) goto 380


      do 140 k=1,nunbond
        i=iun(k)
        j=jun(k)

        e_unbond_drv=0

      if(kunbond(k) == 1) then
        xij=(x(i)-x(j))      
        yij=(y(i)-y(j))      
        zij=(z(i)-z(j))
        rij=sqrt(xij**2+yij**2+zij**2)

        if(i.le.npart1 .and. j.le.npart1) then
          e_unbond_drv=ga1_f1+2*ga2_f1*(gQ_f1-gQ0_f1)
        else if(i.le.npart1 .and. j.gt.npart1) then
          e_unbond_drv = ga1_b + 2 * ga2_b * ( gQ_b - gQ0_b ) + ga1_w + 2 * ga2_w * ( gQ_w - gQ0_w )
        else if(iFlagMov.eq.2) then
          e_unbond_drv=ga1_f2+2*ga2_f2*(gQ_f2-gQ0_f2)
        endif
        
          if((rij/runbond_nat(k)-gr0) .lt. -8*gdr .or. (rij/runbond_nat(k)-gr0) .gt. 8*gdr) then
            e_unbond_drv=0.0
          else
            tem=exp((rij/runbond_nat(k)-gr0)/gdr)
            e_unbond_drv=e_unbond_drv*tem/(1+tem)**2/(runbond_nat(k)*gdr)
          endif

        if(i.le.npart1 .and. j.gt.npart1) then
          e_unbond_drv = e_unbond_drv+alpha_Qw*(ga1_w+2*ga2_w*(gQ_w-gQ0_w))
        endif

        if(i.le.npartM) then
          fxun(i)=fxun(i)+e_unbond_drv*xij/rij
          fyun(i)=fyun(i)+e_unbond_drv*yij/rij
          fzun(i)=fzun(i)+e_unbond_drv*zij/rij
        endif
        if(j.le.npartM) then
          fxun(j)=fxun(j)-e_unbond_drv*xij/rij
          fyun(j)=fyun(j)-e_unbond_drv*yij/rij
          fzun(j)=fzun(j)-e_unbond_drv*zij/rij
        endif
      endif

  140 continue
  380 continue

      return
      end






!  *****************************************************
      SUBROUTINE FUNBOND_with(e_unbond_tot,e_bind_tot)
!  *****************************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)
      common/forceunbonded/fxun(MAXN),fyun(MAXN),fzun(MAXN)
      common/nativeinfo/rbond_nat(MAXN),runbond_nat(MAXUNBO),theta_nat(MAXN),dihedral_nat(MAXN)
      common/constants/pi,boltz,avsn0,gamma,amass,sigma_ij
      common/variables/temp,dt,nstep,nadim,nsnap,gm,enscale
      common/contactmap/iun(MAXUNBO),jun(MAXUNBO),kunbond(MAXUNBO),nQnative_f1,nQnative_f2,nQnative_b
      common/interacparam/ck_r,ck_tht,ck_phi1,ck_phi3,epsil1,epsil2,epsil
      common/LagrangeData/ga1_f1, ga2_f1, ga1_f2, ga2_f2, ga1_b, ga2_b, gQ0_f1, gQ0_f2, gQ0_b, gr0, gdr, gQ_f, gQ_f1, gQ_f2, gQ_b, eGr
      common/LagrangeData2/ga1_w, ga2_w, gQ0_w, gQ_w, alpha_Qw

      common/Solvation/k_sol, n_sol, m_sol, epsilon_p, epsilon_pp, pseudoQ_f,pseudoQ_b, dr_sol, ddr_sol

      common/bias/Alpha1,Alpha2,Beta,Delta
      common/non_native/CritR_non, gQ_non_f, gQ_non_b
      
      
      gQ_f=0.0
      gQ_f1=0.0
      gQ_f2=0.0
      gQ_b=0.0
      gQ_w=0.0
      sum_rij=0.0
      gQ_non_f=0.0
      gQ_non_b=0.0
 

!     r''-r'
      dr_sol=3.0
!     r*-r'
      ddr_sol=1.5


!  unbonded force for A:      
      do 40 k=1,nunbond
        i=iun(k)
        j=jun(k)

          e_unbond=0
          e_unbond_drv=0
      
          xij=x(i)-x(j)      
          yij=y(i)-y(j)      
          zij=z(i)-z(j)
          
          rij=sqrt(xij**2+yij**2+zij**2)

! cancel for Q_w
!          if ( kunbond(k) == 1 .and. rij >= (2*runbond_nat(k)+3.0) .or. &

           if ( kunbond(k) == 2 .and. rij >= (2*CritR_non+3.0     ) .or. & 
               kunbond(k) == 0 .and. rij >= (1.8*sigma_ij        ) .or. &
               i>npartM .and. j>nparM )    goto 60
          
          if ( kunbond(k) == 1 ) then
            if((rij/runbond_nat(k)-gr0) .lt. 8*gdr) then
              tem=1/(1+exp((rij/(runbond_nat(k)+ddr_sol)-gr0)/gdr))
              if(i.le.npart1 .and. j.le.npart1) then
                gQ_f1=gQ_f1+tem
              else if(i.le.npart1 .and. j.gt.npart1) then
                gQ_b=gQ_b+tem
              else if(iFlagMov.eq.2) then
                gQ_f2=gQ_f2+tem
              endif
            endif
            if(i.le.npart1 .and. j.gt.npart1) sum_rij=sum_rij+(rij-runbond_nat(k))

            if( rij .lt. runbond_nat(k) ) then
              Zr=(runbond_nat(k)/rij)**k_sol
              dZr=-k_sol*Zr/rij
              e_unbond     = epsil1*Zr*(Zr-2.0)
              e_unbond_drv = epsil1*2*(Zr-1)*dZr
			             
            else if(rij .lt. runbond_nat(k)+ddr_sol) then
              Yr=(rij-(runbond_nat(k)+ddr_sol))**2
              Yrn=Yr**n_sol
              dYr=2*(rij-(runbond_nat(k)+ddr_sol))
              dr2n=ddr_sol**(2*n_sol)
!  Note that C should be .../(r*-r')^4n
              CC=4*n_sol*(epsil1+epsilon_pp)/dr2n**2
              e_unbond=CC*Yrn*(0.5*Yrn-dr2n)/(2*n_sol)+epsilon_pp
              e_unbond_drv=CC/(2*n_sol)*(Yrn-dr2n)*(n_sol*Yrn/Yr)*dYr

            else            
              Yr=(rij-(runbond_nat(k)+ddr_sol))**2
              dYr=2*(rij-(runbond_nat(k)+ddr_sol))
              ddr=dr_sol-ddr_sol
              BB=epsilon_p*m_sol*ddr**(2*(m_sol-1))
              h1=(1-1.0/m_sol)*ddr**2/(epsilon_p/epsilon_pp+1.0)
              h2=(m_sol-1)*ddr**(2*m_sol)/(1.0+epsilon_pp/epsilon_p)
              e_unbond=-BB*(Yr-h1)/(Yr**m_sol+h2)
              e_unbond_drv=-BB*dYr/(Yr**m_sol+h2)+BB*(Yr-h1)/(Yr**m_sol+h2)**2*m_sol*Yr**(m_sol-1)*dYr
            endif
!  e_unbond_drv=-dU/dr
            e_unbond_drv=-e_unbond_drv

           if(i.le.npart1 .and. j.le.npart1) then
             e_unbond = Alpha1 * e_unbond
             e_unbond_drv = Alpha1 * e_unbond_drv
           else if(i.le.npart1 .and. j.gt.npart1) then
             e_unbond = Beta * e_unbond
             e_unbond_drv = Beta * e_unbond_drv
           else if(iFlagMov.eq.2) then
             e_unbond = Alpha2 * e_unbond
             e_unbond_drv = Alpha2 * e_unbond_drv
           endif     
            
      else if(kunbond(k) == 2) then
         if((rij/CritR_non-gr0) .lt. 8*gdr) then
            tem=1/(1+exp((rij/CritR_non-gr0)/gdr))
            if(i.le.npart1 .and. j.gt.npart1) then
              gQ_non_b=gQ_non_b+tem
            else if(iFlagMov.eq.2 .or. j.le.npart1) then
              gQ_non_f=gQ_non_f+tem
            endif
         endif      	
      	
         rij12=(sigma_ij/rij)**12.0
         rHP  = exp( -1 * (rij - CritR_non )**2 / 2 )  
         if(i.le.npart1 .and. j.le.npart1) then
           Atem = Alpha1
         else if(i.le.npart1 .and. j.gt.npart1) then
!           Atem = sqrt(Alpha1*Alpha2)
           Atem = Beta
         else
           Atem = Alpha2
          endif     
         e_unbond = Atem * epsil2 * rij12 - Delta * epsil1 * rHP
         e_unbond_drv = Atem * 12 * epsil2 * rij12 / rij - Delta * epsil1 * rHP * ( rij - CritR_non )    

!         print *, e_unbond, e_unbond_drv 
      
      else
         rij12=(sigma_ij/rij)**12.0
         if(i.le.npart1 .and. j.le.npart1) then
           Atem = Alpha1
         else if(i.le.npart1 .and. j.gt.npart1) then
!           Atem = sqrt(Alpha1*Alpha2)
           Atem = Beta
         else
           Atem = Alpha2
          endif     
         e_unbond = Atem * epsil2 * rij12
         e_unbond_drv = Atem * 12 * epsil2 * rij12 / rij
      endif

          
          e_unbond_tot=e_unbond_tot+e_unbond
          if(i.le.npart1 .and. j.gt.npart1) e_bind_tot=e_bind_tot+e_unbond
          
        if(i.le.npartM) then
          fxun(i)=fxun(i)+e_unbond_drv*xij/rij
          fyun(i)=fyun(i)+e_unbond_drv*yij/rij
          fzun(i)=fzun(i)+e_unbond_drv*zij/rij
        endif
        if(j.le.npartM) then
          fxun(j)=fxun(j)-e_unbond_drv*xij/rij
          fyun(j)=fyun(j)-e_unbond_drv*yij/rij
          fzun(j)=fzun(j)-e_unbond_drv*zij/rij
        endif
          
   60 continue
   40 continue

      if(iFlagMov.ne.2) gQ_f2=0
      gQ_f=gQ_f1+gQ_f2



!  Lagrange constrant potential: 
!  Fixed B

      gQ_w=gQ_b-alpha_Qw*sum_rij

      eGr_f1 = ga1_f1 * ( gQ_f1 - gQ0_f1 ) + ga2_f1 * ( gQ_f1 - gQ0_f1 ) ** 2
      eGr_f2 = ga1_f2 * ( gQ_f2 - gQ0_f2 ) + ga2_f2 * ( gQ_f2 - gQ0_f2 ) ** 2
      eGr_b = ga1_b * ( gQ_b - gQ0_b ) + ga2_b * ( gQ_b - gQ0_b ) ** 2
      eGr_f=eGr_f1+eGr_f2
      eGr_w = ga1_w * ( gQ_w - gQ0_w ) + ga2_w * ( gQ_w - gQ0_w ) ** 2
      eGr   = eGr_f + eGr_b +eGr_w

      if(eGr.lt.1e-5 .and. eGr.gt.-1e-5) goto 380

      do 340 k=1,nunbond
        i=iun(k)
        j=jun(k)

        e_unbond_drv=0

      if(kunbond(k) == 1) then
        xij=(x(i)-x(j))      
        yij=(y(i)-y(j))      
        zij=(z(i)-z(j))
        rij=sqrt(xij**2+yij**2+zij**2)

        if(i.le.npart1 .and. j.le.npart1) then
          e_unbond_drv=ga1_f1+2*ga2_f1*(gQ_f1-gQ0_f1)
        else if(i.le.npart1 .and. j.gt.npart1) then
          e_unbond_drv = ga1_b + 2 * ga2_b * ( gQ_b - gQ0_b ) + ga1_w + 2 * ga2_w * ( gQ_w - gQ0_w )
        else
          e_unbond_drv=ga1_f2+2*ga2_f2*(gQ_f2-gQ0_f2)
        endif
        
          if((rij/(runbond_nat(k)+ddr_sol)-gr0) .lt. -8*gdr .or. (rij/(runbond_nat(k)+ddr_sol)-gr0) .gt. 8*gdr) then
            e_unbond_drv=0.0
          else
            tem=exp((rij/(runbond_nat(k)+ddr_sol)-gr0)/gdr)
            e_unbond_drv=e_unbond_drv*tem/(1+tem)**2/((runbond_nat(k)+ddr_sol)*gdr)
          endif

        if(i.le.npart1 .and. j.gt.npart1) then
          e_unbond_drv = e_unbond_drv+alpha_Qw*(ga1_w+2*ga2_w*(gQ_w-gQ0_w))
        endif
          
        if(i.le.npartM) then
          fxun(i)=fxun(i)+e_unbond_drv*xij/rij
          fyun(i)=fyun(i)+e_unbond_drv*yij/rij
          fzun(i)=fzun(i)+e_unbond_drv*zij/rij
        endif
        if(j.le.npartM) then
          fxun(j)=fxun(j)-e_unbond_drv*xij/rij
          fyun(j)=fyun(j)-e_unbond_drv*yij/rij
          fzun(j)=fzun(j)-e_unbond_drv*zij/rij
        endif
      endif

  340 continue

  380 continue


      return
      end




!************************************
      SUBROUTINE RANTERM
!************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/randomterms/c_0,c_1,c_2,randconst
      common/forcerandnew/frandx(MAXN),frandy(MAXN),frandz(MAXN)
      do 10 i=1,npartM
        frandx(i)=gauss(xsi)*randconst
        frandy(i)=gauss(xsi)*randconst
        frandz(i)=gauss(xsi)*randconst
  10  continue
      return
      end



!  ************************************************
      SUBROUTINE VERLET(enerkin,e_pot)
!  ***********************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond

	  parameter(nbinmax=105)
      parameter(nEbinmax=105)
	  parameter(nRbinmax=105)

      common/coordinates/x(MAXN),y(MAXN),z(MAXN)
      common/velocity/vx(MAXN),vy(MAXN),vz(MAXN)
   
      common/forcetotalnew/fx(MAXN),fy(MAXN),fz(MAXN)
      common/forcetotalold/fxo(MAXN),fyo(MAXN),fzo(MAXN)
      common/forcerandnew/frandx(MAXN),frandy(MAXN),frandz(MAXN)
      common/forcerandold/frandxo(MAXN),frandyo(MAXN),frandzo(MAXN)
    
      common/gyrationradius/gyr	  

      common/variables/temp,dt,nstep,nadim,nsnap,gm,enscale
      common/randomterms/c_0,c_1,c_2,randconst
      
	  common/interacparam/ck_r,ck_tht,ck_phi1,ck_phi3,epsil1,epsil2,epsil
      common/constants/pi,boltz,avsn0,gamma,amass,sigma_ij
      common/LagrangeData/ga1_f1, ga2_f1, ga1_f2, ga2_f2, ga1_b, ga2_b, gQ0_f1, gQ0_f2, gQ0_b, gr0, gdr, gQ_f, gQ_f1, gQ_f2, gQ_b, eGr
      common/LagrangeData2/ga1_w, ga2_w, gQ0_w, gQ_w, alpha_Qw

      common/Solvation/k_sol, n_sol, m_sol, epsilon_p, epsilon_pp, pseudoQ_f,pseudoQ_b, dr_sol, ddr_sol

!  histogram of (Q_f,Q_b)
      common/HistogramData_fb/nbin_f,nbin_b,dbin_f,dbin_b,vbin0,nbinsnap,nbinsnap0
      common/HistogramData_FB/PFBbin(nbinmax,nbinmax), PFbin(nbinmax), PBbin(nbinmax)

!  histogram of (Q_f,Q_b,E)
      common/EQ_fbHistogramData/nEbin,dEbin,vEbin0, IsEbin
      common/EQ_FBHistogramData/PFBEbin(nbinmax,nbinmax,nEbinmax)

!  histogram of (Q_f,Q_b,E)
      common/EbQbHistogramParameter/nEbbin,dEbbin,vEbbin0, IsEbbin
      common/EbQbHistogramData/PQbEbbin(nbinmax,nEbinmax)

!  histogram of (Q_f,Q_b,R)
      common/RQ_fbHistogramData/nRbin,dRbin,vRbin0, IsRbin
      common/RQ_FBHistogramData/PFBRbin(nbinmax,nbinmax,nRbinmax)

!  histogram of (Q_f1, Q_f2, Q_b)
      common/HistogramData_FFB/PFFBbin(nbinmax,nbinmax,nbinmax)

!  histogram of (Q_w), (Q_w, Q_b) for construction of F(Q_b), (Q_w, R, Qb=0?) for F(R) and K, (Q_w, E_b, Qb=0?) for K changing with epsilon_b
      common/Q_wHistogramParam/nWbin, dWbin, vWbin0, IsWbin, cri_Qb
      common/Q_wHistogramData/PQwbin(nbinmax), PQwQbbin(nbinmax,nbinmax), PQwQfbin(nbinmax,nbinmax), &
                  PQwRbin(nbinmax,nRbinmax,2), PQwEbbin(nbinmax,nEbinmax,2)

	  common/pbc/pL,dl,R

      common/non_native/CritR_non, gQ_non_f, gQ_non_b


      enerkin=0

      call RANTERM

      do 2 j=1,npartM
      x(j)=(x(j)+dt*vx(j)+dt**2*(fxo(j)+frandxo(j)-gm*vx(j))/2.)
      y(j)=(y(j)+dt*vy(j)+dt**2*(fyo(j)+frandyo(j)-gm*vy(j))/2.)
      z(j)=(z(j)+dt*vz(j)+dt**2*(fzo(j)+frandzo(j)-gm*vz(j))/2.)
   2  continue


      call force(e_pot,e_unbond_tot,e_bind_tot,e_tors_tot,e_bend_tot,e_bond_tot)

      do 3 j=1,npartM
      vx(j)=c_1*vx(j)+c_2*(fxo(j)+frandxo(j)+fx(j)+frandx(j))
      vy(j)=c_1*vy(j)+c_2*(fyo(j)+frandyo(j)+fy(j)+frandy(j))
      vz(j)=c_1*vz(j)+c_2*(fzo(j)+frandzo(j)+fz(j)+frandz(j))
      enerkin=enerkin+vx(j)**2+vy(j)**2+vz(j)**2
  3   continue

      do 319 j=1,npartM
      fxo(j)=fx(j) 
      fyo(j)=fy(j) 
      fzo(j)=fz(j) 
      frandxo(j)=frandx(j)
      frandyo(j)=frandy(j)
      frandzo(j)=frandz(j)
  319 continue        
  
      enerkin=0.50*enerkin*amass
      tempins=2.0*enerkin/(3.0*npartM*boltz) 

!  gyr is gyration of radius
!      call calculate_gyrationradius

!  periodic boundary condiction is used based on the coordinate of residue 15
      if(iFlagMov.eq.2 .and. mod(nadim,20).eq.0) call origin_adjust
      if(mod(nadim,20).eq.0) call pbc_shift


!  histogram calculation 
      if(nadim.ge.nbinsnap0 .and. mod(nadim, nbinsnap).eq.0) then

          ibin_f=(gQ_f-vbin0)/dbin_f+1
          if(ibin_f.le.0) ibin_f=1
          if(ibin_f.gt.nbin_f) ibin_f=nbin_f

          ibin_f1=(gQ_f1-vbin0)/dbin_f+1
          if(ibin_f1.le.0) ibin_f1=1
          if(ibin_f1.gt.nbin_f) ibin_f1=nbin_f

          ibin_f2=(gQ_f2-vbin0)/dbin_f+1
          if(ibin_f2.le.0) ibin_f2=1
          if(ibin_f2.gt.nbin_f) ibin_f2=nbin_f

          ibin_b=(gQ_b-vbin0)/dbin_b+1
          if(ibin_b.le.0) ibin_b=1
          if(ibin_b.gt.nbin_b) ibin_b=nbin_b

		  iEbin=(e_pot-vEbin0)/dEbin+1
		  if(iEbin.le.0) iEbin=1
		  if(iEbin.gt.nEbin) iEbin=nEbin

		  iEbbin=(e_bind_tot-vEbbin0)/dEbbin+1
		  if(iEbbin.le.0) iEbbin=1
		  if(iEbbin.gt.nEbbin) iEbbin=nEbbin

		  iRbin=(R-vRbin0)/dRbin+1
		  if(iRbin.le.0) iRbin=1
		  if(iRbin.gt.nRbin) iRbin=nRbin

          iWbin=(gQ_w-vWbin0)/dWbin+1
          if(iWbin.le.0) iWbin=1
          if(iWbin.gt.nWbin) iWbin=nWbin
          if(gQ_b.le.cri_Qb) then
            iQb_cri=1
          else
            iQb_cri=2
          endif

!  for PFBbin:
          PFbin(ibin_f)=PFbin(ibin_f)+1
          PBbin(ibin_b)=PBbin(ibin_b)+1
          PFBbin(ibin_f,ibin_b)=PFBbin(ibin_f,ibin_b)+1
		  PFBEbin(ibin_f,ibin_b,iEbin)=PFBEbin(ibin_f,ibin_b,iEbin)+1
		  PQbEbbin(ibin_b,iEbbin)=PQbEbbin(ibin_b,iEbbin)+1
		  PFBRbin(ibin_f,ibin_b,iRbin)=PFBRbin(ibin_f,ibin_b,iRbin)+1
		  PFFBbin(ibin_f1,ibin_f2,ibin_b)=PFFBbin(ibin_f1,ibin_f2,ibin_b)+1

          PQwbin(iWbin)=PQwbin(iWbin)+1
          PQwQbbin(iWbin,ibin_b)=PQwQbbin(iWbin,ibin_b)+1
          PQwRbin(iWbin,iRbin,iQb_cri)=PQwRbin(iWbin,iRbin,iQb_cri)+1
          PQwEbbin(iWbin,iEbbin,iQb_cri)=PQwEbbin(iWbin,iEbbin,iQb_cri)+1
          PQwQfbin(iWbin,ibin_f1)=PQwQfbin(iWbin,ibin_f1)+1
        
      endif


	  
	  if(nadim.eq.0)          write(9,'(''        nadim      gQ_f      gQ_b       gQ_w       E_k       E_pot    E_b      E_bind       eGr         R'')')
	  if(nadim.eq.0)          write(*,'(''        nadim      gQ_f      gQ_b       gQ_w       E_k       E_pot    E_b      E_bind       eGr         R'')')
	  if(mod(nadim,nsnap)==0) write(9,'(i13,10f10.3)') nadim, gQ_f, gQ_b, gQ_w, enerkin, e_pot, e_bond_tot, e_unbond_tot, e_bind_tot, eGr, R
	  if(mod(nadim,nsnap)==0) write(*,'(i13,10f10.3)') nadim, gQ_f, gQ_b, gQ_w, enerkin, e_pot, e_bond_tot, e_unbond_tot, e_bind_tot, eGr, R

	
      return 
      end



!  ********************************************
       SUBROUTINE PBC_SHIFT
!  *******************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)
      common/pbc/pL,dl,R

	  if(npart1.le.0) goto 400


! note that the chain 2 is located at the center

	  x00=0.0
	  y00=0.0
	  z00=0.0
	  do 10 i=1,npart1
	    x00=x00+x(i)
		y00=y00+y(i)
		z00=z00+z(i)
 10   enddo

      x00=x00/npart1
      y00=y00/npart1
      z00=z00/npart1

!  R is the distance to the center (0,0,0) of the box to center of chain A (x00,y00,z00) 
	  R=sqrt(x00**2+y00**2+z00**2)
      

	  ixflag = 0   
	  if ( x00 > ( pL / 2.0 + dl ) ) then
	    ixflag = 1
	    do 2 i = 1, npart1
	      x(i) = x(i) - pL
2       continue
      endif

	  if ( ixflag == 1 ) goto 111  
	     
	  if ( x00 < ( -1 * pL / 2.0 - dl ) ) then
	    do 12 i = 1, npart1
	      x(i) = x(i) + pL
12      continue
      endif
      
111	  continue

	  iyflag = 0      
	  if ( y00 > ( pL / 2.0 + dl ) ) then
	    iyflag = 1
	    do 22 i = 1, npart1
	      y(i) = y(i) - pL
22      continue
      endif

	  if ( iyflag == 1 ) goto 222
	        
	  if ( y00 < ( -1 * pL / 2.0 - dl ) ) then
	    do 32 i = 1, npart1
	      y(i) = y(i) + pL
32      continue
      endif
      
222   continue

	  izflag = 0      
	  if ( z00 > ( pL / 2.0 + dl ) ) then	  
	    izflag = 1
	    do 42 i = 1, npart1
	      z(i) = z(i) - pL
42      continue
      endif

	  if ( izflag == 1 ) goto 333
	        
	  if ( z00 < ( -1 * pL / 2.0 - dl ) ) then
	    do 52 i = 1, npart1
	      z(i) = z(i) + pL
52      continue
      endif
333   continue	  

400   continue	  
	  	  
	  end





!  ********************************************
      BLOCK DATA
!  *******************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/constants/pi,boltz,avsn0,gamma,amass,sigma_ij
      common/variables/temp,dt,nstep,nadim,nsnap,gm,enscale
      common/gausscoeff/a1,a3,a5,a7,a9

!  PHYSICAL CONSTANTS
!       data pi,boltz,avsn0/3.141592654,1.3806E-23,6.0222E+23/
      data pi/3.141592654/
      data a1/3.949846138/,a3/0.252408784/,a5/0.076542912/,a7/0.008355968/,a9/0.029899776/
      end


!  *********************************************
      SUBROUTINE FTORSION(e_tors_tot)
!  *********************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/bias/Alpha1,Alpha2,Beta,Delta
      parameter(eps=1.0e-3,eps1=1.0-1.0e-6,eps2=1.0e-4) 
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)
      common/forcetorsion/fxph(MAXN),fyph(MAXN),fzph(MAXN)
      common/nativeinfo/rbond_nat(MAXN),runbond_nat(MAXUNBO),theta_nat(MAXN),dihedral_nat(MAXN)
      common/interacparam/ck_r,ck_tht,ck_phi1,ck_phi3,epsil1,epsil2,epsil
      common/constants/pi,boltz,avsn0,gamma,amass,sigma_ij

      do 30 i=1,npartM-3
      j=i+1
      k=i+2
      l=i+3
      xij=(x(i)-x(j))
      yij=(y(i)-y(j))
      zij=(z(i)-z(j))

      xkj=(x(k)-x(j))
      ykj=(y(k)-y(j))
      zkj=(z(k)-z(j))

      xkl=(x(k)-x(l))
      ykl=(y(k)-y(l))
      zkl=(z(k)-z(l))

      rij= sqrt( xij**2 + yij**2 + zij**2 ) 
      rkj= sqrt( xkj**2 + ykj**2 + zkj**2 ) 
      rkl= sqrt( xkl**2 + ykl**2 + zkl**2 ) 


      xmj=yij*zkj-ykj*zij
      ymj=zij*xkj-zkj*xij
      zmj=xij*ykj-xkj*yij

      xnk=ykj*zkl-ykl*zkj
      ynk=zkj*xkl-zkl*xkj
      znk=xkj*ykl-xkl*ykj

      xil=ymj*znk-ynk*zmj
      yil=zmj*xnk-znk*xmj
      zil=xmj*ynk-xnk*ymj

      rnk=sqrt(xnk**2+ynk**2+znk**2)
      rmj=sqrt(xmj**2+ymj**2+zmj**2)

      if(rnk**2.lt.eps.or.rmj**2.lt.eps) go to 30

      phi=(xnk*xmj+ynk*ymj+znk*zmj)/(rnk*rmj)

         if (abs(phi) .gt. eps1) then
  
         x(l)=x(l)+eps2
         y(l)=y(l)+eps2
         z(l)=z(l)+eps2
         
         xkl=(x(k)-x(l))
         ykl=(y(k)-y(l))
         zkl=(z(k)-z(l))

         rkl= sqrt( xkl**2 + ykl**2 + zkl**2 )
         
         xnk=ykj*zkl-ykl*zkj
         ynk=zkj*xkl-zkl*xkj
         znk=xkj*ykl-xkl*ykj
         
         xil=ymj*znk-ynk*zmj
         yil=zmj*xnk-znk*xmj
         zil=xmj*ynk-xnk*ymj
         
         rnk=sqrt(xnk**2+ynk**2+znk**2)
         rmj=sqrt(xmj**2+ymj**2+zmj**2)
     
         if(rnk**2.lt.eps.or.rmj**2.lt.eps) go to 20
         phi=(xnk*xmj+ynk*ymj+znk*zmj)/(rnk*rmj)

         end if

	  if(phi.gt.eps1) phi=eps1
	  if(phi.lt.-1.0*eps1) phi=eps1

      phi=acos(phi)
      phi=sign(phi,xkj*xil+ykj*yil+zkj*zil)
      phi_0=dihedral_nat(i)


      e_tors=ck_phi1*(1.0-cos(phi-phi_0))+ck_phi3*(1.0-cos(3.0*(phi-phi_0)))

      drv1=ck_phi1*sin(phi-phi_0)

      drv2=3.0*ck_phi3*(4.0*cos(phi)**2.0*sin(phi-3.0*phi_0)+3.0*cos(phi)*sin(3.0*phi_0)-sin(phi)*cos(3.0*phi_0))

      e_tors_drv=-(drv1+drv2)
      
      
      if(i.le.npart1-3) then
        tem=Alpha1
      else if(i.gt.npart1) then
        tem=Alpha2
      else
        tem=0
      endif

      e_tors_tot=e_tors_tot+tem*e_tors

      rijrkj=(xij*xkj+yij*ykj+zij*zkj)/rkj**2
      rklrkj=(xkl*xkj+ykl*ykj+zkl*zkj)/rkj**2

      drix= xmj*rkj/rmj**2
      driy= ymj*rkj/rmj**2
      driz= zmj*rkj/rmj**2

      drlx=-xnk*rkj/rnk**2 
      drly=-ynk*rkj/rnk**2 
      drlz=-znk*rkj/rnk**2 

      drjx= (rijrkj-1)*drix - rklrkj*drlx
      drjy= (rijrkj-1)*driy - rklrkj*drly
      drjz= (rijrkj-1)*driz - rklrkj*drlz

      drkx= (rklrkj-1)*drlx - rijrkj*drix
      drky= (rklrkj-1)*drly - rijrkj*driy
      drkz= (rklrkj-1)*drlz - rijrkj*driz

      fxph(i)=fxph(i)+tem*e_tors_drv*drix 
      fyph(i)=fyph(i)+tem*e_tors_drv*driy 
      fzph(i)=fzph(i)+tem*e_tors_drv*driz

      fxph(j)=fxph(j)+tem*e_tors_drv*drjx
      fyph(j)=fyph(j)+tem*e_tors_drv*drjy
      fzph(j)=fzph(j)+tem*e_tors_drv*drjz

      fxph(k)=fxph(k)+tem*e_tors_drv*drkx
      fyph(k)=fyph(k)+tem*e_tors_drv*drky
      fzph(k)=fzph(k)+tem*e_tors_drv*drkz

      fxph(l)=fxph(l)+tem*e_tors_drv*drlx 
      fyph(l)=fyph(l)+tem*e_tors_drv*drly
      fzph(l)=fzph(l)+tem*e_tors_drv*drlz    


   20 continue
   30 continue

      return
      end


!  *****************************************************
      SUBROUTINE FBOND(e_bond_tot)
!  *****************************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      parameter(eps1=1.0e-6) 
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)
      common/forcestretching/fxr(MAXN),fyr(MAXN),fzr(MAXN)
      common/nativeinfo/rbond_nat(MAXN),runbond_nat(MAXUNBO),theta_nat(MAXN),dihedral_nat(MAXN)
      common/interacparam/ck_r,ck_tht,ck_phi1,ck_phi3,epsil1,epsil2,epsil

  	  
	  do 10 i=1,npartM-1
        if(i.eq.npart1) go to 10
      j=i+1
      xij=(x(i)-x(j))
      yij=(y(i)-y(j))
      zij=(z(i)-z(j))

      rij= sqrt( xij**2 + yij**2 + zij**2 ) 

      if(rij.lt.eps1) go to 10
      e_bond=ck_r*(rij-rbond_nat(i))**2
      
      e_bond_tot=e_bond_tot+e_bond
      e_bond_drv=-2*ck_r*(rij-rbond_nat(i))

      drix=xij/rij
      driy=yij/rij
      driz=zij/rij
      drjx=-drix
      drjy=-driy
      drjz=-driz

      if(i.ne.npart1) then
        fxr(i)=fxr(i)+e_bond_drv*drix
        fyr(i)=fyr(i)+e_bond_drv*driy
        fzr(i)=fzr(i)+e_bond_drv*driz
        fxr(j)=fxr(j)+e_bond_drv*drjx
        fyr(j)=fyr(j)+e_bond_drv*drjy
        fzr(j)=fzr(j)+e_bond_drv*drjz
      endif
   10 continue      
      return
      end


!  ***************************************************
      SUBROUTINE FBEND(e_bend_tot)
!  **************************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/bias/Alpha1,Alpha2,Beta,Delta
      parameter(eps=1.0e-3,epstht=1.0-1.0e-12) 
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)
      common/forcebending/fxth(MAXN),fyth(MAXN),fzth(MAXN)
      common/nativeinfo/rbond_nat(MAXN),runbond_nat(MAXUNBO),theta_nat(MAXN),dihedral_nat(MAXN)
      common/interacparam/ck_r,ck_tht,ck_phi1,ck_phi3,epsil1,epsil2,epsil

      do 30 i=1,npartM-2
      j=i+1
      k=i+2
      xij=(x(i)-x(j))
      yij=(y(i)-y(j))
      zij=(z(i)-z(j))

      xkj=(x(k)-x(j))
      ykj=(y(k)-y(j))
      zkj=(z(k)-z(j))
   
      rij= sqrt( xij**2 + yij**2 + zij**2 ) 
      rkj= sqrt( xkj**2 + ykj**2 + zkj**2 ) 

      if(rij.lt.eps.or.rkj.lt.eps) go to 20
      costheta=(xij*xkj+yij*ykj+zij*zkj)/(rij*rkj)
      if(abs(costheta).gt.epstht) costheta=sign(epstht,costheta)

      theta=acos(costheta)

      delta=theta-theta_nat(i)

      if(i.le.npart1-2) then
        tem=Alpha1
      else if(i.gt.npart1) then
        tem=Alpha2
      else
        tem=0
      endif

      e_bend=tem*ck_tht*delta**2
      e_bend_drv=-tem*2*ck_tht*delta

      e_bend_tot=e_bend_tot+e_bend

      sintinv=1.0/sin(theta)

      drix=sintinv*(costheta*xij/rij-xkj/rkj)/rij
      driy=sintinv*(costheta*yij/rij-ykj/rkj)/rij
      driz=sintinv*(costheta*zij/rij-zkj/rkj)/rij
 
      drkx=sintinv*(costheta*xkj/rkj-xij/rij)/rkj 
      drky=sintinv*(costheta*ykj/rkj-yij/rij)/rkj 
      drkz=sintinv*(costheta*zkj/rkj-zij/rij)/rkj


      fxth(i)=fxth(i)+e_bend_drv*drix
      fyth(i)=fyth(i)+e_bend_drv*driy
      fzth(i)=fzth(i)+e_bend_drv*driz

      fxth(j)=fxth(j)+e_bend_drv*(-drix-drkx)
      fyth(j)=fyth(j)+e_bend_drv*(-driy-drky)
      fzth(j)=fzth(j)+e_bend_drv*(-driz-drkz)

      fxth(k)=fxth(k)+e_bend_drv*drkx
      fyth(k)=fyth(k)+e_bend_drv*drky
      fzth(k)=fzth(k)+e_bend_drv*drkz

   20 continue
   30 continue

      return
      end



!  **********************************************
      SUBROUTINE NATIVEINFORMATION
!  **********************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      parameter(eps=1.0e-3,epstht=1.0-1.0e-12) 
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)
      common/nativeinfo/rbond_nat(MAXN),runbond_nat(MAXUNBO),theta_nat(MAXN),dihedral_nat(MAXN)
      common/contactmap/iun(MAXUNBO),jun(MAXUNBO),kunbond(MAXUNBO),nQnative_f1,nQnative_f2,nQnative_b
      common/constants/pi,boltz,avsn0,gamma,amass,sigma_ij


      nQnative_f1=0
      nQnative_f2=0
	  nQnative_b=0

      do 10 k=1,nunbond
        i=iun(k)
        j=jun(k)
        xij=(x(i)-x(j))
        yij=(y(i)-y(j))
        zij=(z(i)-z(j))
        runbond_nat(k)=sqrt(xij**2+yij**2+zij**2)
        if(kunbond(k).eq.1) then
          if(i.le.npart1.and.j.le.npart1) then
            nQnative_f1=nQnative_f1+1
          else if(i.gt.npart1.and.j.gt.npart1) then
            nQnative_f2=nQnative_f2+1
          else
            nQnative_b=nQnative_b+1
          endif
        endif
  10  continue


      do 9 i=1,nparttol-1
      j=i+1
      xij=(x(j)-x(i))
      yij=(y(j)-y(i))
      zij=(z(j)-z(i))
      rbond_nat(i)=sqrt(xij**2+yij**2+zij**2)
    9 continue

      do 11 i=1,nparttol-2
      j=i+1
      k=i+2
      xij=(x(i)-x(j))
      yij=(y(i)-y(j))
      zij=(z(i)-z(j))

      xkj=(x(k)-x(j))
      ykj=(y(k)-y(j))
      zkj=(z(k)-z(j))

      rij= sqrt( xij**2 + yij**2 + zij**2 ) 
      rkj= sqrt( xkj**2 + ykj**2 + zkj**2 ) 

      if(rij.lt.eps.or.rkj.lt.eps) go to 11

      costheta=(xij*xkj+yij*ykj+zij*zkj)/(rij*rkj)
      if(abs(costheta).gt.epstht) costheta=sign(epstht,costheta)

      theta_nat(i)=acos(costheta)
   11 continue

      do 12 i=1,nparttol-3
      j=i+1
      k=i+2
      l=i+3
      xij=(x(i)-x(j))
      yij=(y(i)-y(j))
      zij=(z(i)-z(j))

      xkj=(x(k)-x(j))
      ykj=(y(k)-y(j))
      zkj=(z(k)-z(j))

      xkl=(x(k)-x(l))
      ykl=(y(k)-y(l))
      zkl=(z(k)-z(l))

      rij= sqrt( xij**2 + yij**2 + zij**2 ) 
      rkj= sqrt( xkj**2 + ykj**2 + zkj**2 ) 
      rkl= sqrt( xkl**2 + ykl**2 + zkl**2 ) 

      xmj=yij*zkj-ykj*zij
      ymj=zij*xkj-zkj*xij
      zmj=xij*ykj-xkj*yij

      xnk=ykj*zkl-ykl*zkj
      ynk=zkj*xkl-zkl*xkj
      znk=xkj*ykl-xkl*ykj

      xil=ymj*znk-ynk*zmj
      yil=zmj*xnk-znk*xmj
      zil=xmj*ynk-xnk*ymj

      rnk=sqrt(xnk**2+ynk**2+znk**2)
      rmj=sqrt(xmj**2+ymj**2+zmj**2)

      if(rnk**2.lt.eps.or.rmj**2.lt.eps) go to 12

      phi=(xnk*xmj+ynk*ymj+znk*zmj)/(rnk*rmj)

      phi=acos(phi)
      dihedral_nat(i)=sign(phi,xkj*xil+ykj*yil+zkj*zil)
   12 continue

      return
      end




!  ***************************************************
      FUNCTION RANNYU()
!  ***************************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(ooto12=1.0/4096.0)
      parameter(itwo12=4096)
      common /rnyucm/ m1,m2,m3,m4,l1,l2,l3,l4
      i1=l1*m4+l2*m3+l3*m2+l4*m1
      i2=l2*m4+l3*m3+l4*m2
      i3=l3*m4+l4*m3
      i4=l4*m4
      l4=mod(i4,itwo12)
      i3=i3+i4/itwo12  
      l3=mod(i3,itwo12)
      i2=i2+i3/itwo12
      l2=mod(i2,itwo12)
      l1=mod(i1+i2/itwo12,itwo12)
      rannyu=ooto12*(float(l1)+ooto12*(float(l2)+ooto12*(float(l3)+   ooto12*(float(l4)))))
      return
      end

!  *****************************************************
      SUBROUTINE SETRN(iseed)  
!  ****************************************************
      implicit real*8 (a-h,o-z) 
      implicit integer*4 (i-n)

      common /rnyucm/ m(4),l(4)
      integer iseed(4)
      do 7 j = 1, 4
      isn = 0
      do 5 i = 1, 4
      ipe = 4 - i
      ipd = 10 ** ipe
      id = iseed(j) / ipd   
      isn = isn + id * 8 ** ipe
 5    iseed(j) = iseed(j) - id * ipd
 7    iseed(j) = isn
         do 10 i=1,4
         l(i)=iseed(i)
   10    continue
      l(4)=2*(l(4)/2)+1
      return
      end



!  *****************************************************
      BLOCK DATA RNYUBD
!  ****************************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      common /rnyucm/ m(4),l(4)
      data m / 502,1521,4071,2107/
      data l /   0,   0,   0,   1/
      end   
         


!  *****************************************************
      SUBROUTINE SAVERN(iseed)
!  *****************************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      common /rnyucm/ m(4),l(4) 
      integer iseed(4)
         do 10 i=1,4  
         iseed(i)=l(i)
   10    continue
      do 30 j = 1,4
      isn = 0
      do 20 i = 1,4  
      ipe = 4 - i
      ipo = 8 ** ipe
      id = iseed(j) / ipo
      isn = isn + id * 10 ** ipe
 20   iseed(j) = iseed(j) - ipo * id
 30   iseed(j) = isn  
      return
      end


!************************************
      FUNCTION GAUSS(xsi)
!************************************
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      common/gausscoeff/a1,a3,a5,a7,a9
      r=0
      do 2 j=1,12 
      r=r+rannyu()
  2   continue
      r=(r-6.0)/4.0
      r2=r**2
      gauss=((((a9*r2+a7)*r2+a5)*r2+a3)*r2+a1)*r
      return
      end



!  *******************************************
      SUBROUTINE calculate_gyrationradius
!  *******************************************
      implicit real*8 (a-h,o-z) 
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/coordinates/x(MAXN),y(MAXN),z(MAXN)
      common/gyrationradius/gyr


      ctx=0
      cty=0
      ctz=0
      anr=0

      do 2113 i=1,nparttol
      ctx=ctx+x(i) 
      cty=cty+y(i) 
      ctz=ctz+z(i) 
      anr=anr+x(i)**2+y(i)**2+z(i)**2
 2113 continue

      ctx=ctx/nparttol
      cty=cty/nparttol
      ctz=ctz/nparttol
      gyr1=anr/nparttol

      gyr=sqrt(gyr1-(ctx**2+cty**2+ctz**2))       

      return
      end 


!  *******************************************
      SUBROUTINE InitVel(enerkin)
!  *******************************************
!     initiliaze velocities
! ccccccccccccccccccccccccccccccccc

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter(MAXN=430,MAXUNBO=(MAXN)*(MAXN-1)/2)
      common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
      common/velocity/vx(MAXN),vy(MAXN),vz(MAXN)
      common/constants/pi,boltz,avsn0,gamma,amass,sigma_ij
      common/variables/temp,dt,nstep,nadim,nsnap,gm,enscale
      common/randomterms/c_0,c_1,c_2,randconst
      common/interacparam/ck_r,ck_tht,ck_phi1,ck_phi3, epsil1,epsil2,epsil

      sumvx=0.0 
      sumvy=0.0 
      sumvz=0.0 
      vm=sqrt(temp*boltz/amass)

      do 1001 i=1,npartM
      vx(i)=vm*gauss(xsi)
      vy(i)=vm*gauss(xsi)
      vz(i)=vm*gauss(xsi)
      sumvx=sumvx+vx(i) 
      sumvy=sumvy+vy(i) 
      sumvz=sumvz+vz(i) 
 1001 continue

      sumvx=sumvx/npartM
      sumvy=sumvy/npartM
      sumvz=sumvz/npartM
      enerkin=0.0
 
      do 3 i=1,npartM
      vx(i)=vx(i)-sumvx 
      vy(i)=vy(i)-sumvy
      vz(i)=vz(i)-sumvz
      enerkin=enerkin+vx(i)**2+vy(i)**2+vz(i)**2
   3  continue
 
      enerkin =  0.5 * enerkin * amass 
      tempins=2.0*enerkin/(3.0*npartM*boltz) 

      return
      end

