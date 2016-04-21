c	base4.f
C
C       MULTISCALE MIXED FINITE ELEMENT BASIS in Square Mesh
C       NODE = 4 OVER-SAMPLING 
C       TO SOLVE THE MATRIX PROBLEM, USE dmgd9v
C
C***************************************************
      program main
      parameter (max_nm=5600594, ! nm for level 10
     +           max_nxf=2049,  ! nxf for level 10
     +           max_nyf=2049)  ! nyf for level 10

      real*8 ma(max_nm*9),mrhs(max_nm),a(max_nxf,max_nyf,9),
     +       vb(max_nm),ldu(3*max_nm),work(max_nxf*12),
     +       wa(max_nm),wb(max_nm),p(max_nm),ppp(max_nm),
     +       perm(-511:2560,-511:2560), 
     +       erroru(2,max_nxf,max_nyf),rhs(max_nxf,max_nyf)

      real*8 hx,hy,tol,resno,area,x,y,eu,eflux,value
      real*8 stf1(4,4),stf2(4,4),stf3(4,4),stf4(4,4)
      real*8 xmax,ymax,x0,y0,xblock,yblock
      real*8 coor,basis,exact
      integer ielement,k1,k2
      integer nx,ny,nxc,nyc,nxf,nyf,nm,levels,i,j,k,l,
     +        istart,iprep,maxit,ifail,index,kcoor,node,
     +        iout(6),n1,n2,n3,iteration,num,it,nxyf,mx,
     +        levelsglobal,ineps,mi,mj,mimax,mjmax,nx0,ny0,
     +        nperm1,nperm2,nperm3,ntran
	  real*8 eps
	  common /epsilon/eps
      common /domain/ xmax,ymax,x0,y0
      common /param/ nx,ny,hx,hy
      data nxc,nyc,istart,iprep,maxit,tol,ifail
     +    /5,5,0,0,100,1.0d-8,110/
      data iout/0,0,0,0,0,0/

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

      ! Over-sampling parameters
      kx=nx
      ky=ny
      nx=4*nx
      ny=4*ny
      xblock=xmax
      yblock=ymax
      xmax=4.0d0*xmax
      ymax=4.0d0*ymax
      levels=levels+2

C.....MULTIGRID INFORMATION...................

      nxf=nxc
      nyf=nyc
      nm=nxf*nyf
      do 10 k=2,levels
        nxf=2*nxf-1
        nyf=2*nyf-1
        nm=nm+nxf*nyf
   10 continue

      print*, 'nxf=',nxf,'nyf=',nyf,'nm=',nm
      print*, 'xmax=',xmax,'hx=',hx
c*********************************************
c	input the permeability data
c*********************************************
c	if (mx*kx .eq. 2048) then
c                nxstart=-511
c                nystart=-511
c        else if (mx*kx .eq. 1024) then	
c		nxstart=-255
c		nystart=-255
c	else if (mx*kx .eq. 512 ) then
c		 nxstart=-127
c                 nystart=-127
c        else if (mx*kx .eq. 256 ) then
c                 nxstart=-63
c                 nystart=-63
c        else if (mx*kx .eq. 64) then
c                 nxstart=-31
c                 nystart=-31
c	else 
c		print*, 'permeability data error'
c		stop
c	endif
c
c	open(1,file='../dat/perms.dat')
c	    read(1,*) nx1,ny1,ntot1
cc
cc.	data check
cc	
c	    n1=2*(nxstart-1)
c	    if (mx*kx .ne. (nx1+n1)) then
c		print*, 'gird not consistent'
c		stop
c	    endif
c	    do i=nxstart,nx1+nxstart-1
c        	read(1,*) (perm(i,j),j=nystart,ny1+nystart-1)
c	    enddo
c	close(1)

       open(1,file='base/base4',status='unknown')
C***********************************************************
      do 3000 mj=1,mjmax
      do 3000 mi=1,mimax
C***********************************************************
      
      x0=xleft+(mi-1)*xblock-1.5d0*xblock ! (x0,y0) is the left-bottom
      y0=yleft+(mj-1)*yblock-1.5d0*yblock       ! corner of the sample
c      nx0=(mi-1)*kx-2*kx+1          
c      ny0=(mj-1)*ky-2*ky+1             
c      print*, 'mi,mj,nx0,ny0,x0,y0=',mi,mj,nx0,ny0,x0,y0

c     do i=nx0,nx0+nx-1
c        do j=ny0,ny0+ny-1
c           perm(i,j)=(coef((i-1)*hx,(j-1)*hy)+coef((i-1)*hx,j*hy)
c    +                +coef(i*hx,(j-1)*hy)+coef(i*hx,j*hy))/4.0d0
c        end do
c     end do

C.....FORM THE MA MATRIX.....................

      do j=1,nyf
         do i=1,nxf
            do k=1,9
               a(i,j,k)=0.0d0
            end do
         end do
      end do

      ! boundary nodes: left and right
      do j=1,ny+1
         a(1,j,5)=1.0d0
         a(nx+1,j,5)=1.0d0
      end do

      ! Boundary nodes: top and bottom
      do i=1,nx+1
         a(i,1,5)=1.0d0
         a(i,ny+1,5)=1.0d0
      end do

      ! Interior nodes
      do j=2,ny
         do i=2,nx
            call get_stiff(stf1,i-1,j-1)
            call get_stiff(stf2,i,j-1)
            call get_stiff(stf3,i-1,j)
            call get_stiff(stf4,i,j)
            if(i.eq.2.or.j.eq.2) then
               a(i,j,1)=0.0d0
            else
               a(i,j,1)=stf1(1,3)
            end if
            if(j.eq.2) then
               a(i,j,2)=0.0d0
            else
               a(i,j,2)=stf1(2,3)+stf2(1,4)
            end if
            if(i.eq.nx.or.j.eq.2) then
               a(i,j,3)=0.0d0
            else
               a(i,j,3)=stf2(2,4)
            end if
            if(i.eq.2) then
               a(i,j,4)=0.0d0
            else
               a(i,j,4)=stf1(4,3)+stf3(1,2)
            end if
            a(i,j,5)=stf1(3,3)+stf2(4,4)+stf3(2,2)+stf4(1,1)
            if(i.eq.nx) then
               a(i,j,6)=0.0d0
            else
               a(i,j,6)=stf2(3,4)+stf4(2,1)
            end if
            if(i.eq.2.or.j.eq.ny) then
               a(i,j,7)=0.0d0
            else
               a(i,j,7)=stf3(4,2)
            end if
            if(j.eq.ny) then
              a(i,j,8)=0.0d0
            else
              a(i,j,8)=stf3(3,2)+stf4(4,1)
            endif
            if(i.eq.nx.or.j.eq.ny) then
               a(i,j,9)=0.0d0
            else
               a(i,j,9)=stf4(3,1)
            end if
         end do
      end do

      ! change to the multigrid solver format
      do j=1,nyf  
         do i=1,nxf
            do k=1,9
               call assign_ma(ma,i,j,k,a(i,j,k),nxf,nyf)
            end do
         end do
      end do

C.....ASSEMBLE THE RHS...................................

      do j=1,nyf
         do i=1,nxf
            rhs(i,j)=0.0d0
         end do
      end do

      ! boundary nodes: left and right
      do j=1,ny+1
         rhs(1,j)=exact(x0,y0+(j-1)*hy)
         rhs(nx+1,j)=exact(x0+nx*hx,y0+(j-1)*hy)
      end do

      ! Boundary nodes: top and bottom
      do i=1,nx+1
         rhs(i,1)=exact(x0+(i-1)*hx,y0)
         rhs(i,ny+1)=exact(x0+(i-1)*hx,y0+ny*hy)
      end do

      ! Interior nodes
      do j=2,ny
         do i=2,nx
            if(i.eq.2.and.j.ne.2.and.j.ne.ny) then
               call get_stiff(stf1,i-1,j-1)
               call get_stiff(stf3,i-1,j)
               rhs(i,j)=rhs(i,j)
     +            -stf1(1,3)*rhs(i-1,j-1)
     +            -(stf1(4,3)+stf3(1,2))*rhs(i-1,j)
     +            -stf3(4,2)*rhs(i-1,j+1)
            end if
            if(i.eq.nx.and.j.ne.2.and.j.ne.ny) then
               call get_stiff(stf2,i,j-1)
               call get_stiff(stf4,i,j)
               rhs(i,j)=rhs(i,j)
     +                 -stf2(2,4)*rhs(i+1,j-1)
     +                 -(stf2(3,4)+stf4(2,1))*rhs(i+1,j)
     +                 -stf4(3,1)*rhs(i+1,j+1)
            end if
            if(j.eq.2.and.i.ne.2.and.i.ne.nx) then
               call get_stiff(stf1,i-1,j-1)
               call get_stiff(stf2,i,j-1)
               rhs(i,j)=rhs(i,j)
     +            -stf1(1,3)*rhs(i-1,j-1)
     +            -(stf1(2,3)+stf2(1,4))*rhs(i,j-1)
     +            -stf2(2,4)*rhs(i+1,j-1)
            end if
            if(j.eq.ny.and.i.ne.2.and.i.ne.nx) then
               call get_stiff(stf3,i-1,j)
               call get_stiff(stf4,i,j)
               rhs(i,j)=rhs(i,j)
     +            -stf3(4,2)*rhs(i-1,j+1)
     +            -(stf3(3,2)+stf4(4,1))*rhs(i,j+1)
     +            -stf4(3,1)*rhs(i+1,j+1)
            end if
            if(i.eq.2.and.j.eq.2) then
               call get_stiff(stf1,i-1,j-1)
               call get_stiff(stf2,i,j-1)
               call get_stiff(stf3,i-1,j)
               rhs(i,j)=rhs(i,j)
     +            -stf1(1,3)*rhs(i-1,j-1)
     +            -(stf1(2,3)+stf2(1,4))*rhs(i,j-1)
     +            -stf2(2,4)*rhs(i+1,j-1)
     +            -(stf1(4,3)+stf3(1,2))*rhs(i-1,j)
     +            -stf3(4,2)*rhs(i-1,j+1)
            end if
            if(i.eq.2.and.j.eq.ny) then
               call get_stiff(stf1,i-1,j-1)
               call get_stiff(stf3,i-1,j)
               call get_stiff(stf4,i,j)
               rhs(i,j)=rhs(i,j)
     +            -stf1(1,3)*rhs(i-1,j-1)
     +            -(stf1(4,3)+stf3(1,2))*rhs(i-1,j)
     +            -stf3(4,2)*rhs(i-1,j+1)
     +            -(stf3(3,2)+stf4(4,1))*rhs(i,j+1)
     +            -stf4(3,1)*rhs(i+1,j+1)
            end if              
            if(i.eq.nx.and.j.eq.2) then
               call get_stiff(stf1,i-1,j-1)
               call get_stiff(stf2,i,j-1)
               call get_stiff(stf4,i,j)
               rhs(i,j)=rhs(i,j)
     +            -stf1(1,3)*rhs(i-1,j-1)
     +            -(stf1(2,3)+stf2(1,4))*rhs(i,j-1)
     +            -stf2(2,4)*rhs(i+1,j-1)
     +            -(stf2(3,4)+stf4(2,1))*rhs(i+1,j)
     +            -stf4(3,1)*rhs(i+1,j+1)
            end if
            if(i.eq.nx.and.j.eq.ny) then
               call get_stiff(stf2,i,j-1)
               call get_stiff(stf3,i-1,j)
               call get_stiff(stf4,i,j)
               rhs(i,j)=rhs(i,j)
     +            -stf2(2,4)*rhs(i+1,j-1)
     +            -(stf2(3,4)+stf4(2,1))*rhs(i+1,j)
     +            -stf4(3,1)*rhs(i+1,j+1)
     +            -stf3(4,2)*rhs(i-1,j+1)
     +            -(stf3(3,2)+stf4(4,1))*rhs(i,j+1)
            end if

         end do
      end do

      ! change to the multigrid solver format
      do j=1,nyf  
         do i=1,nxf
            call assign_mrhs(mrhs,i,j,rhs(i,j),nxf,nyf)
         end do
      end do

c      print*, 'go to mgd9v'
      call mgd9v(levels,nxc,nyc,nxf,nyf,nm,iout,istart,iprep,
     +     maxit,tol,mrhs,ma,p,vb,work,ldu,wa,wb,resno,ifail)
c      print*,' ifail: ',ifail,' resno: ',resno
c      print*, 'out of mgd9v'

      k=0
      k1=kx*3/2
      k2=ky*3/2
      do j=k2+1,k2+ky+1
         do i=k1+1,k1+kx+1
            k=k+1
            ppp(k)=p((j-1)*(nx+1)+i)
         end do
      end do

      nxyf=(kx+1)*(ky+1)
c      print*,k,nxyf
      write(1,*) nxyf, (ppp(i),i=1,nxyf)

 3000 continue
      close(1)

      end
C     -------------------------------------------------
C     SUBROUTINE ASSIGN_MA
C     -------------------------------------------------
      subroutine assign_ma(ma,i,j,k,value,nxf,nyf)
      integer i,j,k,nxf,nyf
      real*8 ma(nxf,nyf,9),value

      ma(i,j,k)=value

      return
      end

C     -------------------------------------------------
C     SUBROUTINE ASSIGN_MRHS
C     -------------------------------------------------
      subroutine assign_mrhs(mrhs,i,j,value,nxf,nyf)
      integer i,j,nxf,nyf
      real*8 mrhs(nxf,nyf),value

      mrhs(i,j)=value

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
               !x=coor(node,1,i,j)
               !y=coor(node,2,i,j)
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

C*****************************************************
C        INPUT DATA
C*****************************************************
C     -------------------------------------------
C     EXACT SOLUTION
C     -------------------------------------------
      real*8 function exact(x,y)
      common /domain/ xmax,ymax,x0,y0
      real*8 x,y,xmax,ymax,x0,y0,area
      area=xmax*ymax
      exact=-(x-x0-xmax)*(y-y0)/area
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
