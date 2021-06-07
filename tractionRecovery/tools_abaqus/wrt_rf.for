  
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
	parameter(NUMNODE=100000,NELEMS=100000)
	real*8 vol				!,results(NUMNODE,4)
	real*8 nodes(NUMNODE,3),ddxx(NELEMS)
	real*8   cmath(3,3,3,3,NELEMS),cmathh(3,3,3,3),bss(3,3)
	real*8 	 cmat(3,3,3,3),bs(12,3,3),krn(3,3),I4(3,3,3,3)
	integer conectividades(NELEMS,5)
	END MODULE

c******************************************************************************
      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)

      USE CONEC
      INCLUDE 'ABA_PARAM.INC'
      
C
      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE(ARRAY(1),JRRAY(1,1))	
      INTEGER nn,ii,z1
      REAL*8 U2(NUMNODE),V2(NUMNODE),results(NUMNODE,4)
      REAL*8 W2(NUMNODE)
      CHARACTER(256) JOBDIR
      character*276 filename
      INTEGER im,l
	
      LOVRWRT=1
      nn=1
C
      call POSFIL(KSTEP,KINC,ARRAY,JRCD)
      do K1=1,999999
	call DBFILE(0,ARRAY,JRCD)
	if (JRCD.NE.0) GO TO 110
	KEY=JRRAY(1,2)
C
	if (KEY.EQ.104) then
C	
C	  if(nn.eq.numnode) nn=0
	  U2(nn)=ARRAY(4)
	  V2(nn)=ARRAY(5)
	  W2(nn)=ARRAY(6)
	  nn=nn+1
	ENDif
      enddo
c	STOP
C
  110 CONTINUE
C
c      results=0.0
C
	do I=1,nn-1
	    results(I,1)=I

	    results(I,2)=U2(I)

	    results(I,3)=V2(I)

  	    results(I,4)=W2(I)

	enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c SE ESCRIBE EL ARCHIVO PARA EL GID:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	
      
      call GETOUTDIR(JOBDIR,LENJOBDIR)
      filename=' '
        filename(1:lenjobdir)=jobdir(1:lenjobdir)
        filename(lenjobdir+1:lenjobdir+26)='/reac_for.dat'

        open(UNIT=15,file=filename(1:lenjobdir+26), STATUS='REPLACE')
	
c	if(kstep.eq.1) write(15,10)'Sprout - Prescribed Tip Displacements'
c	write(15,10)''
c	write(15,11)'RESULT "DISPLACEMENTS" (micro-meters), Node, U1,U2,U3'
c	write(15,10)'VALUES'	

	do i=1,nn-1
        write(15,16) i,(results(i,j),j=2,4)
	enddo

c	write(15,10)'END VALUES'

10 	FORMAT(A100)
11 	FORMAT(A100,I4,A100)
16 	FORMAT(I6,1X,3(ES14.4,1X))

	close(15)
	return
      end

