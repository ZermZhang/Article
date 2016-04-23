C
C     ASSEMBLE THE MULTISCALE BASIS: OVER-SAMPLING
C
C***************************************************
      program main
      implicit double precision(a-h,o-z)
      parameter (max_nxf=257, max_nyf=257)  
      parameter (max_nm=max_nxf*max_nyf)  
      dimension  p(4,max_nm),
     +        q(4,max_nm),a(4,4),perm(-511:2560,-511:2560)
      dimension stf(4,4),rhs(4),stf1(4,4)

      common /param/ nx,ny,hx,hy
      common /domain/ xmax,ymax,x0,y0
      common /epsilon/eps


C.....PARAMETERS..............................
      open(1,file='input.dat',status='old')
      read(1,*) levels,nx,mx,levelsglobal,ineps,mimax,mjmax
      read(1,*)
      read(1,*) xleft,yleft,xright,yright
      close(1)

      eps=1.d0/ineps
      ny=nx
      my=mx
      xmax=(xright-xleft)/mx
      ymax=(yright-yleft)/my
      hx=xmax/nx
      hy=ymax/ny
      area=hx*hy      ! area of the element
      print*,nx,ny,mx,my,xleft,yleft,xright,yright,xmax,ymax,hx,hy,eps
c*********************************************
c	input the permeability data
c*********************************************
	if (mx*nx .eq. 2048) then
                nxstart=-511
                nystart=-511
        else if (mx*nx .eq. 1024) then	
		nxstart=-255
		nystart=-255
	else if (mx*nx .eq. 512 ) then
		 nxstart=-127
                 nystart=-127
        else if (mx*nx .eq. 256 ) then
                 nxstart=-63
                 nystart=-63
        else if (mx*nx .eq. 128) then
                 nxstart=-31
                 nystart=-31
        else if (mx*nx .eq. 64) then
                 nxstart=-15
                 nystart=-15
	else 
		print*, 'permeability data error'
		stop
	endif

	open(1,file='../dat/perms.dat')
	    read(1,*) nx1,ny1,ntot1
c
c.	data check
c	
	    n1=2*(nxstart-1)
	    if (mx*nx .ne. (nx1+n1)) then
		print*, 'gird not consistent'
		stop
	    endif
	    do i=nxstart,nx1+nxstart-1
        	read(1,*) (perm(i,j),j=nystart,ny1+nystart-1)
	    enddo
	close(1)

       a(1,1)=6.25d0
       a(1,2)=5.0d0
       a(1,3)=4.0d0
       a(1,4)=5.0d0
       a(2,1)=-3.75d0
       a(2,2)=-2.5d0
       a(2,3)=-2.0d0
       a(2,4)=-3.0d0
       a(3,1)=2.25d0
       a(3,2)=1.5d0
       a(3,3)=1.d0
       a(3,4)=1.5d0
       a(4,1)=-3.75d0
       a(4,2)=-3.d0
       a(4,3)=-2.d0
       a(4,4)=-2.5d0

      open(1,file='base/base1',status='old')
      open(2,file='base/base2',status='old')
      open(3,file='base/base3',status='old')
      open(4,file='base/base4',status='old')     
      open(8,file='base/stf',status='unknown')
      open(9,file='base/rhs',status='unknown')
      open(11,file='base/nbase1',status='unknown')
      open(12,file='base/nbase2',status='unknown')
      open(13,file='base/nbase3',status='unknown')
      open(14,file='base/nbase4',status='unknown')

C*******************************************************
      do 3000 mj=1,mjmax 
      do 3000 mi=1,mimax 
C*******************************************************
C        print*,'mi,mj=',mi,mj

C.....READ THE BASIS INFORMATION
      x0=xleft+(mi-1)*xmax
      y0=yleft+(mj-1)*ymax

      read(1,*) nxyf, (p(1,i),i=1,nxyf)
      read(2,*) nxyf, (p(2,i),i=1,nxyf)
      read(3,*) nxyf, (p(3,i),i=1,nxyf) 
      read(4,*) nxyf, (p(4,i),i=1,nxyf) 
      do k=1,4
        a(k,1)=p(k,1)
        a(k,2)=p(k,nx+1)
        a(k,3)=p(k,nxyf)
        a(k,4)=p(k,nxyf-nx)
      end do

      call inverse(a,4)

      do k=1,4
         do i=1,nxyf
            q(k,i)=0.0d0
            do l=1,4
               q(k,i)=q(k,i)+a(k,l)*p(l,i)
            end do
         end do
      end do
      do k=1,4
         do i=1,nxyf
            p(k,i)=q(k,i)
         end do
      end do

      write(11,*) nxyf, (p(1,i),i=1,nxyf)
      write(12,*) nxyf, (p(2,i),i=1,nxyf)
      write(13,*) nxyf, (p(3,i),i=1,nxyf)
      write(14,*) nxyf, (p(4,i),i=1,nxyf)


      do 100 k=1,4
         do 90 l=1,4
            stf(k,l)=0.0d0
            do 80 i=1,nx
               do 70 j=1,ny
                  call get_stiff(stf1,i,j)
                  call get_num(n1,1,i,j)
                  call get_num(n2,2,i,j)
                  call get_num(n3,3,i,j)
                  call get_num(n4,4,i,j)
                  stf(k,l)=stf(k,l)
     +                     +p(k,n1)*p(l,n1)*stf1(1,1)
     +                     +p(k,n1)*p(l,n2)*stf1(1,2)
     +                     +p(k,n1)*p(l,n3)*stf1(1,3)
     +                     +p(k,n1)*p(l,n4)*stf1(1,4)
     +                     +p(k,n2)*p(l,n1)*stf1(2,1)
     +                     +p(k,n2)*p(l,n2)*stf1(2,2)
     +                     +p(k,n2)*p(l,n3)*stf1(2,3)
     +                     +p(k,n2)*p(l,n4)*stf1(2,4)
     +                     +p(k,n3)*p(l,n1)*stf1(3,1)
     +                     +p(k,n3)*p(l,n2)*stf1(3,2)
     +                     +p(k,n3)*p(l,n3)*stf1(3,3)
     +                     +p(k,n3)*p(l,n4)*stf1(3,4)
     +                     +p(k,n4)*p(l,n1)*stf1(4,1)
     +                     +p(k,n4)*p(l,n2)*stf1(4,2)
     +                     +p(k,n4)*p(l,n3)*stf1(4,3)
     +                     +p(k,n4)*p(l,n4)*stf1(4,4)
70               continue
80            continue
90           continue
100         continue

      write(8,*) ((stf(i,j),j=1,4),i=1,4)


         do 200 k=1,4
            rhs(k)=0.d0
            do 180 i=1,nx
              do 170 j=1,ny
               call get_force(value1,1,i,j)
               call get_force(value2,2,i,j)
               call get_force(value3,3,i,j)
               call get_force(value4,4,i,j)
               call get_num(n1,1,i,j)
               call get_num(n2,2,i,j)
               call get_num(n3,3,i,j)
               call get_num(n4,4,i,j)
               rhs(k)=rhs(k)+p(k,n1)*value1+p(k,n2)*value2
     +                 +p(k,n3)*value3+p(k,n4)*value4
170      continue
180      continue
200      continue
      write(9,*) (rhs(i),i=1,4)

 3000 continue

      close(1)
      close(2)
      close(3)
      close(4)
      close(8)
      close(9)
      close(11)
      close(12)
      close(13)
      close(14)

      stop
      end

C     ------------------------------------------------
C     FUNCTION NTRAN
C     ------------------------------------------------
      integer function ntran(i)
      integer i  
      ntran=i 
      return
      end
C****************************************************
C       ELEMENT INFORMATION
C****************************************************
C     -----------------------------------------------
C     SUBROUTINE GET_STIFF
C     -----------------------------------------------
      subroutine get_stiff(stf,i,j)
      common /domain/ xmax,ymax,x0,y0
      common /param/ nx,ny,hx,hy

      real*8 hx,hy,stf(4,4),area,x,y,bk1,bk2,bl1,bl2,xmax,ymax
      real*8 basis_x,basis_y,coor_gauss,d,x0,y0
      integer nx,ny,i,j,node,k,l

      area=hx*hy

      do k=1,4
         do l=1,4
            stf(k,l)=0.0d0
            do node=1,4
               x=coor_gauss(node,1,i,j)
               y=coor_gauss(node,2,i,j)
               bk1=basis_x(k,x,y,i,j)
               bk2=basis_y(k,x,y,i,j)
               bl1=basis_x(l,x,y,i,j)
               bl2=basis_y(l,x,y,i,j)
               stf(k,l)=stf(k,l)+(area/4.0d0)*
     +              (d(1,1,x,y)*bk1*bl1+d(1,2,x,y)*bk1*bl2
     +              +d(2,1,x,y)*bk2*bl1+d(2,2,x,y)*bk2*bl2)
            end do
         end do
      end do

      return
      end
C     -------------------------------------------------
C     SUBROUTINE GET_FORCE
C     --------------------------------------------------
      ! The integral of the force*basis(node) 
      ! over the element (i,j)
      subroutine get_force(value,node,i,j)
      common /domain/ xmax,ymax,x0,y0
      common /param/ nx,ny,hx,hy
      real*8 hx,hy,value,x,y,area,x0,y0,xmax,ymax
      real*8 force,coor_gauss,basis
      integer nx,ny,index,i,j,node

      area=hx*hy

      value=0.0d0
      do k=1,4
         x=coor_gauss(k,1,i,j)
         y=coor_gauss(k,2,i,j)
         value=value+force(x,y)*basis(node,x,y,i,j)/4.0d0
      end do
      value=value*area

      return
      end
C     -------------------------------------------------
C     SUBROUTINE GET_NUM
C     --------------------------------------------------
      subroutine get_num(num,node,i,j)
      common /param/ nx,ny,hx,hy
      real*8 hx,hy
      integer num,node,i,j,nx,ny


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

C     --------------------------------------------
C     FUNCTION COOR
C     --------------------------------------------
      real*8 function coor(node,kcoor,i,j)
      common /domain/ xmax,ymax,xleft,yleft
      common /param/ nx,ny,hx,hy

      real*8 hx,hy,x0,y0,xleft,yleft,xmax,ymax
      integer nx,ny,index,node,kcoor,i,j

      x0=xleft+(i-1)*hx
      y0=yleft+(j-1)*hy

      if     (node.eq.1.and.kcoor.eq.1) then
        coor=x0+hx
      elseif (node.eq.1.and.kcoor.eq.2) then
        coor=y0+hy/2.0d0
      elseif (node.eq.2.and.kcoor.eq.1) then
        coor=x0+hx/2.0d0
      elseif (node.eq.2.and.kcoor.eq.2) then
        coor=y0+hy
      elseif (node.eq.3.and.kcoor.eq.1) then
        coor=x0
      elseif (node.eq.3.and.kcoor.eq.2) then
        coor=y0+hy/2.0d0
      elseif (node.eq.4.and.kcoor.eq.1) then
        coor=x0+hx/2.0d0
      elseif (node.eq.4.and.kcoor.eq.2) then
        coor=y0
      else
        print*, 'problem here in function coor!'
      endif

      return
      end

C     --------------------------------------------
C     FUNCTION COOR
C     --------------------------------------------
      real*8 function coor_gauss(node,kcoor,i,j)
      common /domain/ xmax,ymax,xleft,yleft
      common /param/ nx,ny,hx,hy

      real*8 hx,hy,x0,y0,tmp,xleft,yleft,xmax,ymax
      integer nx,ny,node,kcoor,i,j

      x0=xleft+(i-1)*hx
      y0=yleft+(j-1)*hy
      tmp=1.0d0/dsqrt(3.0d0)

      if     (node.eq.1.and.kcoor.eq.1) then
        coor_gauss=x0+hx*(1.0d0-tmp)/2.0d0
      elseif (node.eq.1.and.kcoor.eq.2) then
        coor_gauss=y0+hy*(1.0d0-tmp)/2.0d0
      elseif (node.eq.2.and.kcoor.eq.1) then
        coor_gauss=x0+hx*(1.0d0+tmp)/2.0d0
      elseif (node.eq.2.and.kcoor.eq.2) then
        coor_gauss=y0+hy*(1.0d0-tmp)/2.0d0
      elseif (node.eq.3.and.kcoor.eq.1) then
        coor_gauss=x0+hx*(1.0d0+tmp)/2.0d0
      elseif (node.eq.3.and.kcoor.eq.2) then
        coor_gauss=y0+hy*(1.0d0+tmp)/2.0d0
      elseif (node.eq.4.and.kcoor.eq.1) then
        coor_gauss=x0+hx*(1.0d0-tmp)/2.0d0
      elseif (node.eq.4.and.kcoor.eq.2) then
        coor_gauss=y0+hy*(1.0d0+tmp)/2.0d0
      else
        print*, 'problem here in function coor!'
      endif

      return
      end

C     ----------------------------------------------
C     FUNCTION BASIS
C     ----------------------------------------------
      real*8 function basis(node,x,y,i,j)
      common /domain/ xmax,ymax,xleft,yleft
      common /param/ nx,ny,hx,hy
      real*8 hx,hy,x,y,x0,y0,area,xleft,yleft,xmax,ymax
      integer nx,ny,node,kcoor,i,j

      x0=xleft+(i-1)*hx
      y0=yleft+(j-1)*hy
      area=hx*hy

      if     (node.eq.1.) then
        basis=(x-x0-hx)*(y-y0-hy)/area
      elseif (node.eq.2) then
        basis=-(x-x0)*(y-y0-hy)/area
      elseif (node.eq.3) then
        basis=(x-x0)*(y-y0)/area
      elseif (node.eq.4) then
        basis=-(x-x0-hx)*(y-y0)/area
      else
        print*, 'problem here in function basis!'
      endif

      return
      end

C     ----------------------------------------------
C     FUNCTION BASIS_X
C     ----------------------------------------------
      real*8 function basis_x(node,x,y,i,j)
      common /domain/ xmax,ymax,xleft,yleft
      common /param/ nx,ny,hx,hy
      real*8 hx,hy,x,y,x0,y0,area,xleft,yleft,xmax,ymax
      integer nx,ny,index,node,kcoor,i,j

      x0=xleft+(i-1)*hx
      y0=yleft+(j-1)*hy
      area=hx*hy

      if     (node.eq.1.) then
        basis_x=(y-y0-hy)/area
      elseif (node.eq.2) then
        basis_x=-(y-y0-hy)/area
      elseif (node.eq.3) then
        basis_x=(y-y0)/area
      elseif (node.eq.4) then
        basis_x=-(y-y0)/area
      else
        print*, 'problem here in function basis_x!'
      endif


      return
      end

C     ----------------------------------------------
C     FUNCTION BASIS_Y
C     ----------------------------------------------
      real*8 function basis_y(node,x,y,i,j)
      common /domain/ xmax,ymax,xleft,yleft
      common /param/ nx,ny,hx,hy
      real*8 hx,hy,x,y,x0,y0,area,xleft,yleft,xmax,ymax
      integer nx,ny,index,node,kcoor,i,j

      x0=xleft+(i-1)*hx
      y0=yleft+(j-1)*hy
      area=hx*hy

      if     (node.eq.1.) then
        basis_y=(x-x0-hx)/area
      elseif (node.eq.2) then
        basis_y=-(x-x0)/area
      elseif (node.eq.3) then
        basis_y=(x-x0)/area
      elseif (node.eq.4) then
        basis_y=-(x-x0-hx)/area
      else
        print*, 'problem here in function basis!'
      endif

      return
      end
C     -------------------------------------------
C     COEFFICIENT MATRIX
C     -------------------------------------------
      real*8 function d(i,j,x,y)
      real*8 x,y,coef
      integer i,j
      if     (i.eq.1.and.j.eq.1) then
        d=coef(x,y)
      elseif (i.eq.1.and.j.eq.2) then
        d=0.0d0
      elseif (i.eq.2.and.j.eq.1) then
        d=0.0d0
      elseif (i.eq.2.and.j.eq.2) then
        d=coef(x,y)
      else
        print*, 'problem here in function d!'
      endif
      return
      end
C*****************************************************
C        INPUT DATA
C*****************************************************

C     -------------------------------------------
C     SOURCE FUNCTION
C     -------------------------------------------
      real*8 function force(x,y)
      real*8 x,y     
      force=-1.0d0
      end

