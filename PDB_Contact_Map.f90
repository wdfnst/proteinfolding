



!*********************************************

	  PROGRAM PDB_CONTACT_MAP

!                      Yongqi Huang      2009-3-29
!        modified by   Huaiqing Cao     2013-10-12
!           to generate runbond_nat of each contact 

!  the residue number is renumbered for calculation convence
!
!*********************************************

	  write(*,'('' This program read the pdb file to produce the comtact map'')')
	  write(*,'('' Make sure that the disordered chain is at the first place'')')
	  write(*,'('' and the ordered chain is at the second place'')')
	  write(*,'('' Remove all the non-coordinate information'')')
	  write(*,'('' such as "TER,TERMINAL,etc" '')')
	  write(*,'(//)')

	  call RENUMBER_TWO_CHAINS
	  call CONTACT_MAP
 	  
	  end PROGRAM     




!*********************************************

	  SUBROUTINE RENUMBER_TWO_CHAINS

!     renumber residue 
!*********************************************

	  common/CHAINRES/NumResA, NumResB
      parameter ( ResType = 20, MaxAtom = 50000, MaxRes = 500 )
	  character*30 filename, Temp, Atom, Res0, Chain0
      character*3 Res3(ResType), Res(MaxRes)
      character*1 Res1(ResType), Chain(MaxAtom)
      character*1 Chain1, ChainR
      integer NumAtom(MaxAtom), NumRes(MaxAtom)
      integer ResChange

      Res3 = (/"GLY", "ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "ASP", "HIS", "ASN", "GLU", "LYS", "GLN", "MET", "ARG", "SER", "THR", "CYS", "PRO"/)
      Res1 = (/"G", "A", "V", "L", "I", "F", "W", "Y", "D", "H", "N", "E", "K", "Q", "M", "R", "S", "T", "C", "P"/)

	  write ( *, '('' Enter input PDB filename'')')
	  read ( *, * ) filename
	  open ( 10, file = filename, status = 'old' )
	  write ( *, '('' Enter filename to save the renumbered PDB file'')')
	  read ( *, * ) filename
	  open ( 20, file = filename, status = "unknown" )
      write ( *, * ) "Enter filename to save the fasta sequence:"
      read ( *, * ) filename
      open ( 201, file = filename, status = "unknown" )

      read ( 10, 101 ) Temp, NumAtom0, Atom, Res0, Chain0, NumRes0, x, y, z, Beta, Occu
      do k = 1, ResType
        if ( Res0 .eq. Res3(k)) then
            write ( 201, '(a1,$)') Res1(k)
            exit
        endif
        exit
      enddo
      rewind ( 10 )
      
      read ( 10, 101, END = 100 ) Temp, NumAtom0, Atom, Res0, Chain0, NumRes0, x, y, z, Beta, Occu
	  rewind ( 10 )
	  Chain1 = Chain0
	  NumRes1 = 0
	  NumRes2 = NumRes0 - 1

      ResChange = 0
      ChainR = Chain1

	  do i = 1, 10000
	    read ( 10, 101, END = 100 ) Temp, NumAtom(i), Atom, Res(i), Chain(i), NumRes(i), x, y, z, Beta, Occu
!        write ( *,'(i4, 3a4, i4)') NumAtom(i), Atom, Res(i), Chain(i), NumRes(i)
!        if ( i .eq. 50 ) stop
        
        if ( Temp(1:4) .eq. "ATOM" ) then
		    NumAtom(i) = i
            
		    if ( NumRes(i) .ne. NumRes2 ) then
		        NumRes2 = NumRes(i)
                NumRes(i) = NumRes1 + 1
		        NumRes1 = NumRes(i)
                do k = 1, ResType
                    if ( Res(i) .eq. Res3(k)) then
                        write ( 201, '(a1,$)') Res1(k)
                        exit
                    endif
                enddo
		    else
		        NumRes(i) = NumRes1
		    endif
        
            if ( ChainR .ne. Chain(i)) then
                write ( 201, * ) " "
                ChainR = Chain(i)
            endif
        
!        write ( *, '(3a4)') Chain1(1:1), Chain(i), Chain(i)(1:1)
		    if ( Chain1(1:1) .eq. Chain(i)(1:1)) then
		        write ( 20, 101 ) Temp, NumAtom(i), Atom, Res(i), Chain(i), NumRes(i), x, y, z, Beta, Occu
		        NumResA = NumRes(i)
		    else
		        write ( 20, 101 ) Temp, NumAtom(i), Atom, Res(i), Chain(i), NumRes(i) + 20, x, y, z, Beta, Occu
		        NumResB = NumRes(i) + 20
            endif
		endif
    
      enddo
 100  continue
 101  format ( a6, i5, 1x, a4, 1x, a3, 1x, a1, i4, 4x, 3f8.3, 2f6.2 )
	  write ( *, * ) NumResA, NumResB

!      do j = 1, i
!        write ( *, '(i4,a4)') j, Chain(j)
!      enddo
      
!! write the fasta sequence
!      ResChange = 1
!!      ChainR = Chain(1)(1:1)
!      do j = 1, i
!        ResJ = Res(j)(1:3)
!        write ( *, '(a4, i4, a4)') ResJ, ResChange, Chain(j)
!        if ( ResChange .eq. 1 ) then
!            ResR = ResJ
!            do k = 1, ResType
!              write ( *, '(2a4)') ResJ, Res3(k)
!              if ( ResJ .eq. Res3(k) ) then
!                  write ( 201, '(a1,$)') Res1(k)
!                  exit
!              endif
!            enddo
!            ResChange = 0
!        else if ( ResChange .eq. 0 ) then
!            if ( ResJ .ne. ResR ) then
!                ResChange = 1
!            endif
!        endif
!
!        ChainJ = Chain(j)(1:1)
!        write ( *, '(2a3)') ChainJ, ChainR
!        if ( ChainJ .ne. ChainR ) then
!            write ( 201, * ) " "
!            ChainR = ChainJ
!        endif
!      enddo

	  end 



!*********************************************

	  SUBROUTINE CONTACT_MAP

!*********************************************

	  common/CHAINRES/NumResA,NumResB
	  character*20 Temp, Atom, Res, Chain, filename, Atomtype
	  parameter ( MaxAtom = 50000, MaxRes = 500 )
	  dimension NumAtom(MaxAtom), NumRes(MaxAtom), x(MaxAtom), y(MaxAtom), z(MaxAtom), KatomType(MaxAtom), N1(MaxAtom), N2(MaxAtom), Nflag(MaxRes,MaxRes)
	  dimension runbind_nat(MaxRes, MaxRes), CA_Count(MaxAtom)
	  
	  write ( *,'('' Enter file name to save the native contact map'')')
	  read ( *, * ) filename
	  open ( 30, file = filename )
	  write ( *,'('' Enter file name to save the appliable contact map'')')
	  read ( *, * ) filename
	  open ( 40, file = filename )
	  write ( *,'('' Enter file name to save the alpha C atom coordinate'')')
	  read ( *, * ) filename
	  open ( 50, file = filename )

  	  open ( 60, file = 'Temp.dat' )


	  write ( *,'('' Enter critical distance, in the UNIT of A'')')
	  read ( *, * ) Rcrit  ! Critical distance to check contact between any two non-H atom pairs with a separation of 4 residues
	                   ! Normal Rcrit=4.5A  

	  Num = 0
	  rewind ( 20 )
      
	  do i = 1, MaxAtom
        CA_Count(i) = 0
	    read ( 20,'(a6, i5, 1x, a4, 1x, a3, 1x, a1, i4, 4x, 3f8.3, 2f6.2)',END = 100 ) Temp, NumAtom(i), Atom, Res, Chain, NumRes(i), x(i), y(i), z(i), Beta, Occu
        Num = Num + 1
		
        if ( Atom(2:2) == 'H' ) then
		  KatomType(i) = 1
		else
		  KatomType(i) = 0
		endif
		if ( Atom(2:3) == 'CA' ) then
            write ( 50, '(3f10.3)') x(i), y(i), z(i) 
            CA_Count(i) = 1
        endif
	  enddo
100   continue
 
 
! The fourth version of this program has different calculating results of contact residue pairs with
! original data from Yongqi Huang.
! Runbind_nat is the distance between C_alpha atoms of two contact residues.
!   according to Kaya & Chan, JMB_2003_326

      do i = 1, MaxRes
          do j = 1, MaxRes
              runbind_nat(i,j) = 100.0
          enddo
      enddo
      

      do i = 1, Num
	    if ( KatomType(i) .ne. 1 ) then
	      do j = i + 1, Num	  
            ni = NumRes(i)
            nj = NumRes(j)
		    if ( (nj-ni) .ge. 4 .and. KatomType(j) .eq. 0 ) then
		      R = sqrt ( ( x(i) - x(j))**2 + ( y(i) - y(j))**2 + ( z(i) - z(j))**2 )
              if ( CA_Count(i) == 1 .and. CA_Count(j) == 1 ) then
                  runbind_nat(ni,nj) = R
              endif
              if ( R .le. Rcrit ) then
                  write ( 60, * ) ni,nj
              endif
	        endif
		  enddo
		endif
      enddo		  		 

	  rewind ( 60 )

	  NMax = 0
	  i = 0
	  do 
	    i = i + 1 
		read ( 60, *, END = 200 ) N1(i), N2(i)
		NMax = NMax + 1
	  enddo
 200  continue
 
	  do i = 1, NMax
	    Iflag = 0
		do j = 1, i - 1
 		  if ( N2(j) .eq. N2(i) .and. N1(j) .eq. N1(i) ) then
            Iflag = 1
			exit
		  endif
		enddo
		if ( Iflag .eq. 0 ) then
		  write ( 30, * ) N1(i), N2(i)
		  Nflag(N1(i),N2(i)) = 1
		endif
      enddo

      do k = 1, NumResA
	    do l = k + 4, NumResA
          if ( Nflag(k,l) .eq. 1 ) then
		  	write ( 40, * ) k, l, 1, runbind_nat(k,l)
		  else
		    write ( 40, * ) k, l ,0, runbind_nat(k,l)
		  endif
		enddo
		
		do n = NumResA + 21, NumResB		
          if ( Nflag(k,n) .eq. 1 ) then
		  	write ( 40, * ) k, n, 1, runbind_nat(k,n)
		  else
		    write ( 40, * ) k, n, 0, runbind_nat(k,n)
		  endif
        enddo 
      enddo

      do k = NumResA + 21, NumResB
        do l = k + 4, NumResB
          if ( Nflag(k,l) .eq. 1 ) then
              write ( 40, * ) k, l, 1, runbind_nat(k,l)
          else
              write ( 40, * ) k, l, 0, runbind_nat(k,l)
          endif
        enddo
      enddo

      close ( 60, status = 'delete' ) 
      
      
      end



