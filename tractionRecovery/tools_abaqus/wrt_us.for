  
C --------------------------------------------------- ELEMENTO DE USUARIO --
C
C
C     HOMOGENEIZACION ASINTOTICA 3D BASADA EN IMPLEMENTACION DE J. F. RODRIGUEZ
C
C
C
C  Esta subrutina nos permite definir las siguientes variables:
C  RHS()
C  AMATRX()
C  SVARS ()
C  ENERGY()
C 
C  Variables que pueden ser actualizadas
C  PNEWDT
C
	MODULE CONEC
	parameter(NUMNODE=100000)
	END MODULE

c******************************************************************************
      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)

      USE CONEC
      INCLUDE 'ABA_PARAM.INC'
      
C
      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE(ARRAY(1),JRRAY(1,1))	
      INTEGER nn,ii,z1,mm,ll
      REAL*8 stress(NUMNODE,6),stressp(NUMNODE,3),strainp(NUMNODE,3)
      CHARACTER(256) JOBDIR
      character*276 filename
      INTEGER im,l
	
      LOVRWRT=1
      nn=1
      mm=1
      ll=1
   
      
    
      stress = 0.d0
      stressp = 0.d0
      strainp = 0.d0
C
      call POSFIL(KSTEP,KINC,ARRAY,JRCD)
      do K1=1,999999
	call DBFILE(0,ARRAY,JRCD)
	if (JRCD.NE.0) GO TO 110
	KEY=JRRAY(1,2)
C
      
      if (KEY.EQ.11) then
C	

	  stress(mm,1)=ARRAY(3)
	  stress(mm,2)=ARRAY(4)
	  stress(mm,3)=ARRAY(5)
        stress(mm,4)=ARRAY(6)
        stress(mm,5)=ARRAY(7)
        stress(mm,6)=ARRAY(8)
	  mm=mm+1
      ENDif


      if (KEY.EQ.401) then
C	

	  stressp(ll,1)=ARRAY(3)
	  stressp(ll,2)=ARRAY(4)
	  stressp(ll,3)=ARRAY(5)
	  ll=ll+1
      ENDif

      if (KEY.EQ.403) then
C	

	  strainp(nn,1)=ARRAY(3)
	  strainp(nn,2)=ARRAY(4)
	  strainp(nn,3)=ARRAY(5)
	  nn=nn+1
      ENDif
      
      enddo

C
  110 CONTINUE
C

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c SE ESCRIBE EL ARCHIVO PARA EL GID:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	
      
      call GETOUTDIR(JOBDIR,LENJOBDIR)
      filename=' '
        filename(1:lenjobdir)=jobdir(1:lenjobdir)
 
      
      filename(lenjobdir+1:lenjobdir+26)='/stresses.dat'

      if (kstep.eq.1) then  
        open(UNIT=15,file=filename(1:lenjobdir+26), STATUS='REPLACE')
      endif
      if (kstep.gt.1) then
	  open(UNIT=15,file=filename(1:lenjobdir+22), ACCESS='APPEND')
      endif
        

	do i=1,mm-1
        write(15,17) i,(stress(i,j),j=1,6)
      enddo
      
      close(15)
      
      
      filename(lenjobdir+1:lenjobdir+26)='/stressesp.dat'

      if (kstep.eq.1) then  
        open(UNIT=15,file=filename(1:lenjobdir+26), STATUS='REPLACE')
      endif
      if (kstep.gt.1) then
	  open(UNIT=15,file=filename(1:lenjobdir+22), ACCESS='APPEND')
      endif
        

	do i=1,ll-1
        write(15,16) i,(stressp(i,j),j=1,3)
      enddo
      
      close(15)
      
      
      filename(lenjobdir+1:lenjobdir+26)='/strainsp.dat'

      if (kstep.eq.1) then  
        open(UNIT=15,file=filename(1:lenjobdir+26), STATUS='REPLACE')
      endif
      if (kstep.gt.1) then
	  open(UNIT=15,file=filename(1:lenjobdir+22), ACCESS='APPEND')
      endif
        

	do i=1,nn-1
        write(15,16) i,(strainp(i,j),j=1,3)
      enddo
      
      

10 	FORMAT(A100)
11 	FORMAT(A100,I4,A100)
16 	FORMAT(I6,1X,3(ES14.4,1X))
17 	FORMAT(I6,1X,6(ES14.4,1X))

	close(15)
	return
      end

