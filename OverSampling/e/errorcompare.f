C   *
C                                                  *
C***************************************************
      program main
      implicit double precision(a-h,o-z)
      parameter (max_nx=512, max_ny=512)
      parameter (max_nxfine=2048, max_nyfine=2048)  
      parameter (max_nxf=max_nx+1, max_nyf=max_ny+1)  

      allocatable :: p(:,:,:,:),utmp(:)
      dimension u(max_nxf*max_nyf),
     +          u_ex((max_nxfine+1)*(max_nyfine+1))
      dimension uhf(4),uhm(4),uec(4),bk1(4),bk2(4)

      character*35 filename

      common /param1/nx,ny,hx,hy,xleft,yleft
      common /param2/nxl,nyl,hxl,hyl
C.....PARAMETERS..............................
      open(1, file='input.dat',status='old')
      read(1,*) levelslocal,nxlocal,nxglobal,levelsglobal,ineps,
     +          mimax,mjmax
	read(1,*)
	read(1,*) xleft,yleft,xright,yright
      close(1)
      pi=dacos(-1.d0)

      levels=levelsglobal     ! multigrid levels
      nx=nxglobal        ! number of spaces in x direction
      ny=nx        ! number of spaces in y direction

      hx=(xright-xleft)/nx
      hy=(yright-yleft)/ny

C.....LOCAL PARAMETERS.........................

      nxl=nxlocal
      nyl=nxl
      hxl=hx/nxl
      hyl=hy/nyl
C      area=0.5d0*hxl*hyl   ! area of the fine element
      open(2,file='result/Omsfem_solution.dat')
      read(2,*) nxf,nyf,(u(k),k=1,nxf*nyf)
      close(2)

c	compare with the exact solution
c
      open(1,file='../dat/result/exact_2048.dat',form='formatted')
      read(1,*)nxfine,nyfine,(u_ex(k),k=1,nxfine*nyfine)   
      close(1)
      nxll=(nxfine-1)/nx/nxl
      nyll=(nyfine-1)/ny/nyl
      hxf=hxl/nxll
      hyf=hyl/nyll
      area=hxf*hyf   ! area of the finest grid
C
       nz=(nxl+1)*(nyl+1)
       allocate(p(4,nx,ny,nz))
       allocate(utmp(nz))
C       print*,nx,ny,nz
      open(1,file='base/nbase1',status='old')
      open(2,file='base/nbase2',status='old')
      open(3,file='base/nbase3',status='old')
      open(4,file='base/nbase4',status='old')

C*******************************************************
      do 1000 mj=1,mjmax 
      do 1000 mi=1,mimax 
C*******************************************************
c        print*,'mi,mj=',mi,mj

C.....READ THE BASIS INFORMATION

      read(1,*) nxyf, (p(1,mi,mj,i),i=1,nxyf)
      read(2,*) nxyf, (p(2,mi,mj,i),i=1,nxyf)
      read(3,*) nxyf, (p(3,mi,mj,i),i=1,nxyf) 
      read(4,*) nxyf, (p(4,mi,mj,i),i=1,nxyf)

 1000 continue

      close(1)
      close(2)
      close(3)
      close(4)

      eu=0.d0
      value=0.d0
      enorml2=0.d0
      p_max=0.0d0
      h1error=0.d0
      h1exact=0.d0

c      open(3,file='result/Ovfine_solution.dat')
      do 3000 mj=1,mjmax
      do 3000 mi=1,mimax
        call get_numg(n1,1,mi,mj)
        call get_numg(n2,2,mi,mj)
        call get_numg(n3,3,mi,mj)
        call get_numg(n4,4,mi,mj)

        do j=1,nyl+1
         do i=1,nxl+1
           call get_num(k1,1,i,j)
           utmp(k1)=u(n1)*p(1,mi,mj,k1)+u(n2)
     +        *p(2,mi,mj,k1)+u(n3)*p(3,mi,mj,k1)
     +        +u(n4)*p(4,mi,mj,k1) 
c           write(3,*)utmp(k1)
         end do
        end do
        do 180 i=1,nxl
         ni=(mi-1)*nxl+i
         do 170 j=1,nyl
          nj=(mj-1)*nyl+j
            do node=1,4
              call get_num(n1,node,i,j)
              uhm(node)=utmp(n1)
            enddo
            do 1801 il=1,nxll     
            do 1701 jl=1,nyll
            do nodel=4,1,-1
             call get_numf(ilf,jlf,nodel,mi,mj,i,j,
     &              il,jl,nxll,nyll)
             xf=xleft+(ilf-1)*hxf
             yf=yleft+(jlf-1)*hyf
             uhtmp=0.d0
             do m=1,4
              uhtmp=uhtmp+uhm(m)*basisf(m,xf,yf,ni,nj)! finer mesh corresponding to 1/mx*nx
             enddo
C             print*,'uhtmp=',uhtmp,nxll,nyll
             uhf(nodel)=uhtmp
             uec(nodel)=u_ex((jlf-1)*nxfine+ilf)        
             dctmp=dabs(uhf(nodel)-uec(nodel))
             if(value.lt.dctmp) value=dctmp
             if(p_max.lt.dabs(uec(nodel))) p_max=dabs(uec(nodel))
            enddo
            call get_integral(el2tmp,eutmp,h1tmp,euxtmp,uhf,uec,
     &         hxf,hyf)
            eu=eu+el2tmp
            enorml2=enorml2+eutmp
            h1error=h1error+h1tmp
            h1exact=h1exact+euxtmp
1701       continue
1801      continue
170      continue
180     continue

3000     continue

        enorml2=dsqrt(enorml2)
 	eu=dsqrt(eu)
        h1error=dsqrt(h1error)
        h1exact=dsqrt(h1exact)
cccccccccccccccccccccccccc
c  This creates a file called error????.dat, where ???? is the partition
c  number on x-direction.
cccccccccccccccccccccccccc
        iter=nx
            ii1 = iter / 1000
            m1 = mod(ii1,10)
            ii2 = iter / 100
            m2 = mod(ii2,10)
            ii3 = iter / 10
            m3 = mod(ii3,10)
            m4 = mod(iter,10)

       filename = '../dat/result/Omsfem_error'//char(m1+48)//
     $          char(48+m2)//char(48+m3)//char(48+m4)//'.dat'


        open(12,file=filename)
 	write(12,*) 'ineps is  ', ineps
 	write(12,*) 'coarse mesh size  ', nxglobal
 	write(12,*) 'fine mesh size ', nxlocal
 	write(12,*) 'L^2 error is  ', eu,eu/enorml2
 	write(12,*) "L^\\infty  error is  ",  value,value/p_max
 	write(12,*) "H^1  error is  ",  h1error,h1error/h1exact
        write(12,*) 'l2, L_inf and H1-norm norm of exact solution: ',
     &       enorml2, p_max ,h1exact
 	close(12)

      stop
      end


C     ------------------------------------------------
C     FUNCTION NTRAN
C     ------------------------------------------------
      integer function ntran(i)
      integer i
      real*8  a

      ntran=i

      return
      end
C     -------------------------------------------------
C     SUBROUTINE GET_NUM
C     --------------------------------------------------
      subroutine get_num(num,node,i,j)
      implicit double precision(a-h,o-z)
      common /param2/nx,ny,hx,hy
      
       if     (node.eq.1) then
        num=(j-1)*(nx+1)+i
      elseif (node.eq.2) then
        num=(j-1)*(nx+1)+i+1
      elseif (node.eq.3) then
        num=j*(nx+1)+i+1
      elseif (node.eq.4) then
        num=j*(nx+1)+i
      else
        print*, 'problem here in get_num!'
      endif

      return
      end
C     -------------------------------------------------
C     SUBROUTINE GET_NUM Global
C     --------------------------------------------------
      subroutine get_numg(num,node,i,j)
      implicit double precision(a-h,o-z)
      common /param1/nx,ny,hx,hy,xleft,yleft
      if     (node.eq.1) then
        num=(j-1)*(nx+1)+i
      elseif (node.eq.2) then
        num=(j-1)*(nx+1)+i+1
      elseif (node.eq.3) then
        num=j*(nx+1)+i+1
      elseif (node.eq.4) then
        num=j*(nx+1)+i
      else
        print*, 'problem here in get_numg!'
      endif

      return
      end
C     ----------------------------------------------
C     FUNCTION BASIS of local element
C     ----------------------------------------------
      real*8 function basisf(node,x,y,i,j)
      implicit double precision(a-h,o-z)
      common /param1/nx,ny,hx,hy,xleft,yleft
      common /param2/nxl,nyl,hxl,hyl

      x0=xleft+(i-1)*hxl
      y0=yleft+(j-1)*hyl
      area=hxl*hyl

      if     (node.eq.1.) then
        basisf=(x-x0-hxl)*(y-y0-hyl)/area
      elseif (node.eq.2) then
        basisf=-(x-x0)*(y-y0-hyl)/area
      elseif (node.eq.3) then
        basisf=(x-x0)*(y-y0)/area
      elseif (node.eq.4) then
        basisf=-(x-x0-hxl)*(y-y0)/area
      else
        print*, 'problem here in local function basis!'
      endif
C        print*,'basisf=', basisf
      return
      end

C     ----------------------------------------------
C     FUNCTION BASIS_X
C     ----------------------------------------------
      real*8 function basisff_x(node,x,y,i,j,hxf,hyf)
      implicit double precision(a-h,o-z)
      common /param1/nx,ny,hx,hy,xleft,yleft
      common /param2/nxl,nyl,hxl,hyl

      x0=xleft+(i-1)*hxf
      y0=yleft+(j-1)*hyf
      area=hxf*hyf

      if     (node.eq.1.) then
        basisff_x=(y-y0-hyf)/area
      elseif (node.eq.2) then
        basisff_x=-(y-y0-hyf)/area
      elseif (node.eq.3) then
        basisff_x=(y-y0)/area
      elseif (node.eq.4) then
        basisff_x=-(y-y0)/area
      else
        print*, 'problem here in function basisff_x!'
      endif


      return
      end

C     ----------------------------------------------
C     FUNCTION BASIS_Y
C     ----------------------------------------------
      real*8 function basisff_y(node,x,y,i,j,hxf,hyf)
      implicit double precision(a-h,o-z)
      common /param1/nx,ny,hx,hy,xleft,yleft
      common /param2/nxl,nyl,hxl,hyl

      x0=xleft+(i-1)*hxf
      y0=yleft+(j-1)*hyf
      area=hxf*hyf

      if     (node.eq.1.) then
        basisff_y=(x-x0-hxf)/area
      elseif (node.eq.2) then
        basisff_y=-(x-x0)/area
      elseif (node.eq.3) then
        basisff_y=(x-x0)/area
      elseif (node.eq.4) then
        basisff_y=-(x-x0-hxf)/area
      else
        print*, 'problem here in function basisff_y!'
      endif

      return
      end
C     -------------------------------------------------
C     SUBROUTINE integral on the finest grid using the rule of three middle points
C     --------------------------------------------------
      subroutine get_integral(el2tmp,eutmp,h1tmp,euxtmp,uhf,uec,
     &           hxf,hyf)

      implicit double precision(a-h,o-z)
      dimension uhf(4),uec(4),tmp(4),A(4,4)

      el2tmp=0.d0
      eutmp=0.d0
      h1tmp=0.d0
      euxtmp=0.d0

      area=hxf*hyf

      do m=1,4
        tmp(m)=uhf(m)-uec(m)
      enddo

      do i=1,4
        A(i,i)=1.d0/9.d0
      enddo
      do i=1,3
        A(i,i+1)=1.d0/18.d0
      enddo
      do i=1,2
        A(i,i+2)=1.d0/36.d0
      enddo
      A(1,4)=1.d0/18.d0
      do j=1,4
       do i=j+1,4
         A(i,j)=A(j,i)
       end do
      end do     
 
      do j=1,4
      do i=1,4
       el2tmp=el2tmp+tmp(i)*tmp(j)*A(i,j)
       eutmp=eutmp+uec(i)*uec(j)*A(i,j)
      end do
      end do
      el2tmp=el2tmp*area
      eutmp=eutmp*area

      do i=1,4
        A(i,i)=2.d0/3.d0
      enddo
      do i=1,3
        A(i,i+1)=-1.d0/6.d0
      enddo
      do i=1,2
        A(i,i+2)=-1.d0/3.d0
      enddo
      A(1,4)=-1.d0/6.d0
      do j=1,4
       do i=j+1,4
         A(i,j)=A(j,i)
       end do
      end do     
 
      do j=1,4
      do i=1,4
       h1tmp=h1tmp+tmp(i)*tmp(j)*A(i,j)
       euxtmp=euxtmp+uec(i)*uec(j)*A(i,j)
      end do
      end do
      end
C     -------------------------------------------------
C     SUBROUTINE GET_NUM
C     --------------------------------------------------
      subroutine get_numf(ni,nj,node,mi,mj,i1,j1,
     &                il,jl,nxll,nyll)
      implicit double precision(a-h,o-z)
      common /param1/nx,ny,hx,hy,xleft,yleft
      common /param2/nxl,nyl,hxl,hyl


      i=((mi-1)*nxl+i1-1)*nxll+il
      j=((mj-1)*nyl+j1-1)*nyll+jl  !Global fine mesh location

      if     (node.eq.1) then
        ni=i
        nj=j
      elseif (node.eq.2) then
        ni=i+1
        nj=j
      elseif (node.eq.3) then
        ni=i+1
        nj=j+1
      elseif (node.eq.4) then
        ni=i
        nj=j+1
      else
        print*, 'problem here in get_numf!'
      endif

      return
      end
