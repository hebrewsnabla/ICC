C******************************************************C
C                                                      C
C   TOTAL MO will be decomposed to fragment MO.        C
C            CODED BY IKEDA ATSUSHI  '06 2/13          C
C    LAST MODIFIED BY IKEDA ATSUSHI                    C
C                                                      C
C******************************************************C

C******************************************************C
C                                                      C
C       Decide MAXNUM                                  C
C                                                      C
C******************************************************C

      PROGRAM MAIN

      Character*10 CC,CD,CE,CW,CA,CB
      
      CALL Getarg(1,CC)
      READ(CC,*) MAXNUM

      CALL Getarg(2,CD)
      READ(CD,*) NUMA

      CALL Getarg(3,CE)
      READ(CE,*) NUMB

      CALL Getarg(4,CW)
      READ(CW,*) NUMOCC

      
      CALL MODECP(MAXNUM,NUMA,NUMB,NUMOCC)

      END
      
      SUBROUTINE MODECP(MAXNUM,NUMA,NUMB,NUMOCC)
            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
C**********************************************************************C
C    VARIABLE NAME                                                     C
C       NUMTOT=number of total basis                                   C
C       NUMA  =number of fragment A basis                              C
C       NUMB  =number of fragment B basis                              C
C     CAF(I,J)=expansion coef of TOTAL MOI AOJ                         C
C     CBF(I,J)=expansion coef of FRAGMENT A's MOI AOJ (I,J<NUMA)       C
C                                FRAGMENT B's MOI AOJ (I,J>NUMB)       C
C              0.0 (other)                                             C
C     CWT(I,J)=expansion coef of TOTAL MOI FRAGMENT MOJ                C
C     SAO(I,J)=overlap between AOI and AOJ                             C
C     SMO(I,J)=overlap between TOTAL MOI and fragment MOJ              C
C     SMOX(I,J)=overlap between fragment MOI and fragment MOJ          C
C     POP(I)=Mullikene-like population of fragment MOI                 C
C**********************************************************************C
C    NAMING                                                            C
C      FRAGMENT MOI (I<NUMA -->fragmentA)                              C
C      FRAGMENT MOI (I>NUMB -->fregmentB+NUMA)                         C
C**********************************************************************C

      DIMENSION CAF(MAXNUM,MAXNUM),CBF(MAXNUM,MAXNUM)
      DIMENSION CWT(MAXNUM,MAXNUM),SMOSUB(MAXNUM,MAXNUM)
      DIMENSION SAO(MAXNUM,MAXNUM),SMO(MAXNUM,MAXNUM)
      DIMENSION SMOX(MAXNUM,MAXNUM)
      DIMENSION IPIV(MAXNUM),CWTSUM(MAXNUM),POP(MAXNUM)
      DIMENSION OCCTOT(MAXNUM), OCC(MAXNUM)
      DIMENSION SMOFRA(MAXNUM,MAXNUM)
      DIMENSION PS(MAXNUM,MAXNUM),P(MAXNUM,MAXNUM)
      DIMENSION OWN(MAXNUM), POL(MAXNUM), CT(MAXNUM)
      DIMENSION OUTV(MAXNUM,4)
C      DIMENSION SCC(MAXNUM,MAXNUM,MAXNUM), SCCSUM(MAXNUM,MAXNUM)
      DIMENSION SCC(MAXNUM,MAXNUM)
      
      CHARACTER*255 FILE1,FILE2,FILE3,FILE4,POINT,KEYWORD,CHAR
      DATA ZERO,ONE,TWO/0.0D+00,1.0D+00,2.0D+00/
!
      Dimension Temp(MAXNUM, MAXNUM)   
C                                                     C
C     ---- Setting the input and output files ----    C
C                                                     C

      NFT1=1
      NFT2=2
      NFT3=3
      NFT4=4

      FILE1='TOTAL'
      FILE2='FRAGA'
      FILE3='FRAGB'
      FILE4='LOG'
      
      OPEN(NFT1,FILE=FILE1)
      OPEN(UNIT=NFT4,FILE=FILE4)

C                                            C
C   ---- input initial information ----      C
C                                            C
       
      NUMTOT=MAXNUM

C                               C
C     ---- Clear Array ----     C
C                               C
      
       DO I=1,MAXNUM
         DO J=1,MAXNUM

       CAF(I,J)=ZERO  
       CBF(I,J)=ZERO
       CWT(I,J)=ZERO
       SAO(I,J)=ZERO
       SMO(I,J)=ZERO
       SMOX(I,J)=ZERO
       
       IF(I.eq.J) THEN
       CWT(I,J)=ONE
       ENDIF
       
         ENDDO

       CWTSUM(I)=ZERO
       POP(I)=ZERO  
       ENDDO

C                                         C
C     ---- READ OVERLAP MARIX ----        C
C                                         C
      
      POINT=' *** Overlap ***'
  601 READ(NFT1,'(255a)') KEYWORD    
      IF(KEYWORD(1:16).eq.POINT(1:16)) GOTO 901
      GOTO 601
      
  901 CONTINUE

      NMAX=NUMTOT/5

      DO N=0,NMAX
       NNOW=N*5+1
       NLAST=NUMTOT-NNOW

       IF(NLAST.eq.-1) GOTO 10

      READ(NFT1,*) CHAR

      READ(NFT1,*) M,SAO(NNOW,NNOW)
      IF(NLAST.eq.0) GOTO 10

      READ(NFT1,*) M,SAO(NNOW+1,NNOW),SAO(NNOW+1,NNOW+1)
      IF(NLAST.eq.1) GOTO 10

      READ(NFT1,*) M,SAO(NNOW+2,NNOW),SAO(NNOW+2,NNOW+1),
     & SAO(NNOW+2,NNOW+2)
      IF(NLAST.eq.2) GOTO 10

      READ(NFT1,*) M,SAO(NNOW+3,NNOW),SAO(NNOW+3,NNOW+1),
     & SAO(NNOW+3,NNOW+2),SAO(NNOW+3,NNOW+3)
      IF(NLAST.eq.3) GOTO 10

      READ(NFT1,*) M,SAO(NNOW+4,NNOW),SAO(NNOW+4,NNOW+1),
     & SAO(NNOW+4,NNOW+2),SAO(NNOW+4,NNOW+3),SAO(NNOW+4,NNOW+4)
      IF(NLAST.eq.4) GOTO 10

      READ(NFT1,*) (M,SAO(I,NNOW),SAO(I,NNOW+1),
     & SAO(I,NNOW+2),SAO(I,NNOW+3),SAO(I,NNOW+4),I=NNOW+5,NUMTOT)

      ENDDO

   10 CONTINUE 

      SUM=ZERO
      
      DO I=1,NUMTOT
        SUM=SUM+SAO(I,2)
      ENDDO

      SUM=ZERO
      

C                                      C
C     ---- READ MO COEFFICIENTS ----   C
C          FROM TOTAL                  C
      
      POINT='     Molecular Orbital Coefficients'
  612 READ(NFT1,"(255a)") KEYWORD    
      IF(KEYWORD(1:35).eq.POINT(1:35)) GOTO 912
      GOTO 612
      
  912 CONTINUE

      NMAX=NUMTOT/5

      DO N=0,NMAX-1
       NNOW=N*5+1

      READ(NFT1,*) CHAR
      READ(NFT1,*) CHAR
      READ(NFT1,*) CHAR
      
      READ(NFT1,1000) (CAF(I,NNOW),CAF(I,NNOW+1),
     & CAF(I,NNOW+2),CAF(I,NNOW+3),CAF(I,NNOW+4),I=1,NUMTOT)

      ENDDO

      READ(NFT1,*) CHAR
      READ(NFT1,*) CHAR
      READ(NFT1,*) CHAR

      NLAST=NUMTOT-NMAX*5
      N=NMAX*5+1     
      
      IF(NLAST.eq.1) THEN 
      READ(NFT1,1001) (CAF(I,N),I=1,NUMTOT)      
      ENDIF
      
      IF(NLAST.eq.2) THEN 
      READ(NFT1,1002) (CAF(I,N),CAF(I,N+1),I=1,NUMTOT)      
      ENDIF
      
      IF(NLAST.eq.3) THEN 
      READ(NFT1,1003) (CAF(I,N),CAF(I,N+1),CAF(I,N+2),I=1,NUMTOT)      
      ENDIF
      
      IF(NLAST.eq.4) THEN 
      READ(NFT1,1004) (CAF(I,N),CAF(I,N+1),CAF(I,N+2),
     & CAF(I,N+3),I=1,NUMTOT)      
      ENDIF
      
       Close(NFT1) 

      
C                                         C
C     ---- READ MO COEFFICIENTS ----      C
C          FROM FRAG A                    C
      
      OPEN(NFT2,FILE=FILE2)

      POINT='     Molecular Orbital Coefficients'
  613 READ(NFT2,"(255a)") KEYWORD    
      IF(KEYWORD(1:35).eq.POINT(1:35)) GOTO 913
      GOTO 613
      
  913 CONTINUE

      NMAX=NUMA/5

      DO N=0,NMAX-1
       NNOW=N*5+1

      READ(NFT2,*) CHAR
      READ(NFT2,*) CHAR
      READ(NFT2,*) CHAR
      
      READ(NFT2,1000) (CBF(I,NNOW),CBF(I,NNOW+1),
     & CBF(I,NNOW+2),CBF(I,NNOW+3),CBF(I,NNOW+4),I=1,NUMA)

      ENDDO

      READ(NFT2,*) CHAR
      READ(NFT2,*) CHAR
      READ(NFT2,*) CHAR

      NLAST=NUMA-NMAX*5
      N=NMAX*5+1     
      
      IF(NLAST.eq.1) THEN 
      READ(NFT2,1001) (CBF(I,N),I=1,NUMA)      
      ENDIF
      
      IF(NLAST.eq.2) THEN 
      READ(NFT2,1002) (CBF(I,N),CBF(I,N+1),I=1,NUMA)      
      ENDIF
      
      IF(NLAST.eq.3) THEN 
      READ(NFT2,1003) (CBF(I,N),CBF(I,N+1),CBF(I,N+2),I=1,NUMA)      
      ENDIF
      
      IF(NLAST.eq.4) THEN 
      READ(NFT2,1004) (CBF(I,N),CBF(I,N+1),CBF(I,N+2),
     & CBF(I,N+3),I=1,NUMA)      
      ENDIF

       Close(NFT2) 

      
C                                         C
C     ---- READ MO COEFFICIENTS ----      C
C          FROM FRAG B                    C
      
      OPEN(NFT3,FILE=FILE3)

      POINT='     Molecular Orbital Coefficients'
  614 READ(NFT3,"(255a)") KEYWORD    
      IF(KEYWORD(1:35).eq.POINT(1:35)) GOTO 914
      GOTO 614
      
  914 CONTINUE

      NMAX=NUMB/5

      DO N=0,NMAX-1
       NNOW=NUMA+N*5+1

      READ(NFT3,*) CHAR
      READ(NFT3,*) CHAR
      READ(NFT3,*) CHAR
      
      READ(NFT3,1000) (CBF(I,NNOW),CBF(I,NNOW+1),
     & CBF(I,NNOW+2),CBF(I,NNOW+3),CBF(I,NNOW+4),I=NUMA+1,NUMTOT)

      ENDDO

      READ(NFT3,*) CHAR
      READ(NFT3,*) CHAR
      READ(NFT3,*) CHAR

      NLAST=NUMB-NMAX*5
      N=NUMA+NMAX*5+1     
      
      IF(NLAST.eq.1) THEN 
      READ(NFT3,1001) (CBF(I,N),I=NUMA+1,NUMTOT)      
      ENDIF
      
      IF(NLAST.eq.2) THEN 
      READ(NFT3,1002) (CBF(I,N),CBF(I,N+1),I=NUMA+1,NUMTOT)      
      ENDIF
      
      IF(NLAST.eq.3) THEN 
      READ(NFT3,1003) (CBF(I,N),CBF(I,N+1),CBF(I,N+2),I=NUMA+1,NUMTOT)
      ENDIF
      
      IF(NLAST.eq.4) THEN 
      READ(NFT3,1004) (CBF(I,N),CBF(I,N+1),CBF(I,N+2),
     & CBF(I,N+3),I=NUMA+1,NUMTOT)      
      ENDIF

       Close(NFT3) 

      Write(NFT4,*) 'read part end'
C                                            C
C    ---- Fill the occupation vector ----    C
C                                            C
 
      DO I=1,NUMTOT
       OCCTOT(I)=ZERO
      ENDDO
      DO I=1,NUMOCC
       OCCTOT(I)=TWO
      ENDDO

      Write(NFT4,*) 'Fill the occupation vector part end'
       
C                                               C
C     ---- Create Full Overlap Matrix ----      C
C                                               C
      
      DO I=1,NUMTOT
        DO J=1,NUMTOT

        IF(I.gt.J) THEN
         SAO(J,I)=SAO(I,J)
        ENDIF

        ENDDO
      ENDDO
      
      Write(NFT4,*) 'Create Full Overlap Matrix part end'
C                                                  C
C     ---- Calculate overlap between MO-MO ----    C
C                                                  C
        
      Write(NFT4,*) 'Convert part start'

!       DO I=1,NUMTOT
!         DO J=1,NUMTOT
!       
!        SUM=ZERO
!       
!         DO K=1,NUMTOT
!             DO L=1,NUMTOT
!             
!        SUM=SUM+CAF(K,I)*CBF(L,J)*SAO(K,L)
!             
!             ENDDO
!           ENDDO  
!
!        SMO(I,J)=SUM
!        SMOSUB(J,I)=SUM
!
!         ENDDO
!       ENDDO

        Do I=1,NUMTOT
           Do J=1,NUMTOT
              Temp(J, I)  = 0.0D0
              SMO(J,I)    = 0.0D0
              SMOSUB(J,I) = 0.0D0
           End do
        End do

       
!       Do I=1,NUMTOT
!          Do J=1,NUMTOT
!             Do L=1,NUMTOT
!                Temp(L,I) = Temp(L,I) + SAO(L,J)*CBF(J,I)
!             End do 
!          End do 
!       End do
!
!       Do I=1,NUMTOT
!          Do J=1,NUMTOT
!             Do L=1,NUMTOT
!                SMO(L,I) = SMO(L,I) + CAF(J,L)*Temp(J,I)
!             End do 
!          End do 
!       End do
!
        N = NUMTOT
        Alpha = 1.0D0
        Beta = 0.0D0
        Call DGemm('N','N',N,N,N,Alpha,SAO,N,CBF,N,Beta,Temp,N)
        Call DGemm('T','N',N,N,N,Alpha,CAF,N,Temp,N,Beta,SMO,N)
        Do I=1,NUMTOT
           Do J=1,NUMTOT 
              SMOSUB(I,J)= SMO(J,I)
           End do
        End do

      Write(NFT4,*) 'Calculate overlap between MO-MO part start'

C                                                    C
C     ---- Calculate overlap between MO-MO FRAGMENT  C
C                                                    C

     
!       DO I=1,NUMTOT
!         DO J=1,NUMTOT
!       
!        SUM=ZERO
!       
!         DO K=1,NUMTOT
!             DO L=1,NUMTOT
!             
!        SUM=SUM+CBF(K,I)*CBF(L,J)*SAO(K,L)
!             
!             ENDDO
!           ENDDO  
!
!        SMOFRA(I,J)=SUM
!
!         ENDDO
!       ENDDO
!     
!       SUM=ZERO
        Call DGemm('T','N',N,N,N,Alpha,CBF,N,Temp,N,Beta,SMOFRA,N)

      Write(NFT4,*) 'Calculate overlap between MO-MO FRAGMENT part end'
C                                          C
C     ---- Solve the equation ----         C
C                                          C
      Write(NFT4,*) 'Enter into LAPACK part'

         CALL DGESV(MAXNUM,MAXNUM,SMO,MAXNUM,IPIV
     &  ,CWT,MAXNUM,INFO)
      Write(NFT4,*) 'LAPACK end'
     
C                                              C
C    ---- Print out the Information ----       C
C                                              C


      WRITE(NFT4,*) INFO
      WRITE(NFT4,'(1x,A30)') "End of Calculation......"

      WRITE(NFT4,'(1x,26A)') "Number of Basis Functions"      
      WRITE(NFT4,'(2x,A5,2x,A5,2x,A5)') "TOTAL","FRAGA","FRAGB"
      WRITE(NFT4,'(2x,I5,2x,I5,2x,I5)') NUMTOT,NUMA,NUMB

C                     
C      MAKING P-MATRIX
C                     
      
!     SUM=ZERO
!     SUMM=ZERO
!     DO K=1,MAXNUM
!     DO L=1,MAXNUM
!       DO I=1,MAXNUM
!         DO J=1,MAXNUM
!           SUM=SUM+CWT(K,I)*CWT(J,I)*SMOFRA(J,L)*OCCTOT(I)
!           SUMM=SUMM+CWT(K,I)*CWT(J,I)*SMOFRA(J,L)*OCCTOT(I)
!         ENDDO                  
!         SCC(K,L,I)=SUMM
!         SUMM=ZERO
!       ENDDO
!      PS(K,L)=SUM
!      SUM=ZERO 
!     ENDDO
!     ENDDO  

      Do I=1,NUMTOT
         Do J=1,NUMTOT
            Temp(J, I)  = 0.0D0
         End do
      End do
      Call DGemm('T','N',N,N,N,Alpha,CWT,N,SMOFRA,N,Beta,Temp,N)

      Do L=1,NUMTOT
         Do I=1,NUMTOT
C            Temp(L, L) = Temp(L, I)*OCCTOT(I)
            Temp(I, L) = Temp(I, L)*OCCTOT(I)
         End do
      End do


C      Do I = 1, NUMTOT
C         Do L = 1, NUMTOT
C            Do K = 1, NUMTOT
C               SCC(K,L,I)=CWT(K,I)*Temp(I,L)
C            End do
C         End do
C      End do

      Do I=1,NUMTOT
         Do J=1,NUMTOT
            SCC(J, I)    = 0.0D0
         End do
      End do

      Do I = 1, NUMTOT
         Do L = 1, NUMTOT
               SCC(L,I)=CWT(L,I)*Temp(I,L)
         End do
      End do
C
      Call DGemm('N','N',N,N,N,Alpha,CWT,N,Temp,N,Beta,PS,N)

C                     
C      MAKING P-MATRIX
C                     
      
      SUM=ZERO
      DO K=1,MAXNUM
      DO L=1,MAXNUM
        DO I=1,MAXNUM
        SUM=SUM+CWT(K,I)*CWT(L,I)*OCCTOT(I)
        ENDDO
      P(K,L)=SUM
      SUM=ZERO
      ENDDO
      ENDDO

      
C                                C
C          POPULATION DECOMPOSE  C
C                                C
      
      DO I=1,MAXNUM
       IF(I.le.NUMA) THEN
        DO J=1,NUMA
         SUM=SUM+PS(I,J)*PS(J,I)/TWO
        ENDDO
        OWN(I)=SUM
        SUM=ZERO

        DO J=NUMA+1,MAXNUM
         SUM=SUM+PS(I,J)*PS(J,I)/TWO
        ENDDO
        CT(I)=SUM
        SUM=ZERO

        ELSE
                
        DO J=1,NUMA
         SUM=SUM+PS(I,J)*PS(J,I)/TWO
        ENDDO
        CT(I)=SUM
        SUM=ZERO

        DO J=NUMA+1,MAXNUM
         SUM=SUM+PS(I,J)*PS(J,I)/TWO
        ENDDO
        OWN(I)=SUM
        SUM=ZERO

        ENDIF
      ENDDO

      DO I=1,MAXNUM
       POP(I)=PS(I,I)
      ENDDO
      
      DO I=1,MAXNUM
       OUTV(I,1)=POP(I)
       OUTV(I,2)=OWN(I)
       OUTV(I,3)=CT(I)
       OUTV(I,4)=OWN(I)+CT(I)
      ENDDO

!      WRITE(NFT4,*) "P-Matrix"
!      DO I=1,MAXNUM
!      DO J=1,MAXNUM
!        QQ=ABS(P(I,J))
!        IF (QQ.gt.0.2) THEN
!           IF (I.le.NUMA) THEN
!             IF (J.le.NUMA) THEN
!           WRITE(NFT4,'(A1,1x,I3,1x,A3,1x,A1,1x,I3,1x,1F10.5)')
!     & "A",I ,"and", "A",J, P(I,J)           
!             ELSE
!          WRITE(NFT4,'(A1,1x,I3,1x,A3,1x,A1,1x,I3,1x,1F10.5)')
!    & "A",I ,"and", "B",J-NUMA, P(I,J)           
!            ENDIF
!           ELSE          
!             IF (J.le.NUMA) THEN
!           WRITE(NFT4,'(A1,1x,I3,1x,A3,1x,A1,1x,I3,1x,1F10.5)')
!     & "B",I-NUMA ,"and", "A",J, P(I,J)           
!             ELSE
!           WRITE(NFT4,'(A1,1x,I3,1x,A3,1x,A1,1x,I3,1x,1F10.5)')
!     & "B",I-NUMA ,"and", "B",J-NUMA, P(I,J)  
!             ENDIF
!           ENDIF 
!         ENDIF  
!      ENDDO
!      ENDDO
      
!      WRITE(NFT4,*) "Orbital Composition"
!      
!      DO I=1,NUMOCC
!        WRITE(NFT4,*) "TOTAL MO", I
!        DO J=1,MAXNUM
!          QT=ABS(CWT(J,I))
!           IF (QT.gt.0.05) THEN
!            IF (J.le.NUMA) THEN
!           WRITE(NFT4,'(A6,2x,A1,I3,2x,1F10.5)')
!     & "FRAGMO","A",J,CWT(J,I)
!           ELSE
!           K=J-NUMA 
!           WRITE(NFT4,'(A6,2x,A1,I3,2x,1F10.5)')
!     & "FRAGMO","B",K,CWT(J,I)
!           ENDIF
!          ENDIF
!        ENDDO
!      ENDDO        
      
      SUM=ZERO
      
      DO I=1,MAXNUM
        DO J=1,MAXNUM
          SUM=SUM+CWT(I,J)*OCCTOT(J)
        ENDDO
       CWTSUM(I)=SUM
      SUM=ZERO
      ENDDO 

      WRITE(NFT4,*) "  "
      WRITE(NFT4,*) "---- Mullikene-like population ----"
      call PRVECZ(OUTV,NUMTOT,NFT4,MAXNUM,NUMA)

      WRITE(NFT4,*) " "
      WRITE(NFT4,*) "Population of FragA"
      
      SUM=ZERO
      Do I=1,NUMA
       SUM=SUM+POP(I)
      ENDDO 
      
      WRITE(NFT4,'(1F10.3)') SUM
      
      WRITE(NFT4,*) "Population of FragB"

      SUM=ZERO
      Do I=NUMA+1,NUMTOT
       SUM=SUM+POP(I)
      ENDDO 
      
      WRITE(NFT4,'(1F10.3)') SUM

      WRITE(NFT4,*) "Population of Total"

      SUM=ZERO
      Do I=1,NUMTOT
       SUM=SUM+POP(I)
      ENDDO 
      
      WRITE(NFT4,'(1F10.3)') SUM

      WRITE(NFT4,*) " "
      WRITE(NFT4,*) "TOTAL MO POPULATION "
      
      DO I=1,NUMOCC
      WRITE(NFT4,*) "TOTAL MO", I
      
        DO K=1,NUMA
         WRITE(NFT4,'(A6,I4,1F10.3)') "FRAG A", K, SCC(K,I)
        ENDDO

        DO L=NUMA+1,NUMTOT
         LL=L-NUMA
         WRITE(NFT4,'(A6,I4,1F10.3)') "FRAG B", LL, SCC(L,I)
        ENDDO

        SUM=ZERO
        DO KK=1,NUMTOT
         SUM=SUM+SCC(KK,I)
        ENDDO
        WRITE(NFT4,'(A9,1F10.3)') "Summation", SUM
        SUM=ZERO
        WRITE(NFT4,*) " "
         
      ENDDO
      
      
      Close(NFT4) 

      RETURN

 1000 FORMAT(21x,5F10.3)
 1001 FORMAT(21x,1F10.3)
 1002 FORMAT(21x,2F10.3)
 1003 FORMAT(21x,3F10.3)
 1004 FORMAT(21x,4F10.3)

      END 

C                                 C
C                                 C
C   PRINT OUT SUBROUTINE          C
C                                 C
C                                 C

      Subroutine PRSQMT(A,NUM,NT,MAXNUM)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MAXNUM,MAXNUM)

      MAX = 5 
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. NUM) IMAX = NUM
      WRITE (NT,9008)
      WRITE (NT,9028) (I,I = IMIN,IMAX)
      WRITE (NT,9008)

      DO 120 J = 1,NUM
  120 WRITE (NT,9048) J,(A(J,I),I = IMIN,IMAX)
      IF (IMAX .LT. NUM) GO TO 100

      RETURN

 9008 FORMAT(1X)
 9028 FORMAT(15X,10(4X,I4,3X))
 9048 FORMAT(I5,2X,8x,10F11.6)

      END

      Subroutine PRSQMZ(A,NUM,NT,MAXNUM,NUMA)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MAXNUM,MAXNUM)
      CHARACTER MOLE

      WRITE (NT,9008)
      WRITE (NT,*) "   TOTAL COMPLEX MO INDEX "
      
      MAX = 5 
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. NUM) IMAX = NUM
      WRITE (NT,9008)
      WRITE (NT,9028) (I,I = IMIN,IMAX)
      WRITE (NT,9008)

      DO 120 J = 1,NUM

      IF(J.LE.NUMA) THEN
      MOLE='A' 
      NUMIN=J
      ELSE
      MOLE='B'
      NUMIN=J-NUMA
      ENDIF
      
  120 WRITE (NT,9048) J,MOLE,NUMIN,(A(J,I),I = IMIN,IMAX)
      IF (IMAX .LT. NUM) GO TO 100

      RETURN

 9008 FORMAT(1x)     
 9028 FORMAT(15X,10(4X,I4,3X))
 9048 FORMAT(I5,2X,1A,1x,I4,2x,10F11.3)

      END
      
      Subroutine PRVECZ(A,NUM,NT,MAXNUM,NUMA)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MAXNUM,4)
      CHARACTER MOLE

      WRITE (NT,9008)
      WRITE (NT,9009) "POP","OWN", "CT", "OWN+CT"
      
      DO 120 J = 1,NUM

      IF(J.LE.NUMA) THEN
      MOLE='A' 
      NUMIN=J
      ELSE
      MOLE='B'
      NUMIN=J-NUMA
      ENDIF
      
  120 WRITE (NT,9048) J,MOLE,NUMIN,A(J,1),A(J,2),A(J,3),A(J,4)

      RETURN

 9008 FORMAT(1x)     
 9009 FORMAT(15X,A11,A11,A11,A11)
 9048 FORMAT(I5,2X,1A,1x,I4,2x,4F11.3)

      END
      
