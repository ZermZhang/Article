c--------------------------------------------------
c***************************************************
C                                                  *
C       SOLVE THE POISSON EQUATION USING THE       *
C       TRIANGLULAR LINEAR FINITE ELEMENT          * 
C       BOUNDARY CONDITION: Dirichlet              *
C       TO SOLVE THE MATRIX PROBLEM, USE dmgd9v    *
C                                                  *
C***************************************************
c	to solve the problem on the  fine grid
c	and  yield the fine grid solution.

c
      program main
      parameter (max_nm=5600594, ! nm for level 9
     +           max_nxf=2049,  ! nxf for level 9
     +           max_nyf=2049)  ! nyf for level 9
  
      real*8 ma(max_nm*9),mrhs(max_nm),a(max_nxf,max_nyf,9),
     +       vb(max_nm),ldu(3*max_nm),work(max_nxf*12),
     +       wa(max_nm),wb(max_nm),p(max_nm),rhs(max_nxf,max_nyf),
     +		perm(-511:2560,-511:2560),
     +       v_1(0:max_nxf-1,max_nyf-1),v_2(max_nxf-1,0:max_nyf-1)
      real*8 hx,hy,tol,resno,area,x,y,eflux,eps,coef
      real*8 stf1(3,3),stf2(3,3),stf3(3,3),stf4(3,3),stf5(3,3),stf6(3,3)
      real*8 value1,value2,value3,value4,value5,value6
      real*8 coor,basis,exact,basis_x,basis_y,r0,delta,pi
      real*8 xleft,xright,yleft,yright,x0,y0,value_x,value_y
      integer nx,ny,nxc,nyc,nxf,nyf,nm,levels,i,j,k,l,
     +        istart,iprep,maxit,ifail,index,kcoor,node,
     +        iout(6),n1,n2,n3,num,ineps
	real*8 d1,d2,d3,d4,value_xp,value_yp,eu,value,norml2,p_max

	character*12 filename
        common/epsilon/eps
      common /param/ nx,ny,hx,hy,xleft,yleft,xright,yright
      data nxc,nyc,istart,iprep,maxit,tol,ifail
     +    /5,5,0,0,50,1.0d-9,110/
      data iout/0,0,0,0,0,0/
	
C.....PARAMETERS..............................
!	open the initial data file 'input.dat'
       open (10,file='../comm/input_lin2048.dat',status='Old')
        read(10,*)
!       print*, 'input nx,ny,levels'
        read(10,*) nx,ny,levels
       close(10)
       open (11,file='../comm/domain.dat',status='Old')
        read(11,*) xleft,yleft,xright,yright
       close(11)

       hx=(xright-xleft)/nx
       hy=(yright-yleft)/ny
       area=0.5d0*hx*hy      ! area of the element
       pi=dacos(-1.0d0)
       open (12,file='../comm/input_eps.dat',status='Old')
        read(12,*) ineps
       close(12)
       eps=1.d0/ineps

C.....MULTIGRID INFORMATION...................

      nxf=nxc
      nyf=nyc
      nm=nxf*nyf
      do k=2,levels
         nxf=2*nxf-1
         nyf=2*nyf-1
         nm=nm+nxf*nyf
      end do

      print*, 'nxf=',nxf,'nyf=',nyf,'nm=',nm
      print*, 'hx=',hx, 'hy=',hy


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
	    value=perm(i-1,j-1)
            call get_stiff(stf1,2,i-1,j-1,value)
	    value=perm(i-1,j-1)
            call get_stiff(stf2,1,i-1,j-1,value)
	    value=perm(i,j-1)
            call get_stiff(stf3,2,i,j-1,value)
	    value=perm(i-1,j)
            call get_stiff(stf4,1,i-1,j,value)
	    value=perm(i,j)
            call get_stiff(stf5,2,i,j,value)
	    value=perm(i,j)
            call get_stiff(stf6,1,i,j,value)
c            if(i.eq.2.and.j.eq.2) then
c               print*, 'stf1='
c               print*, ((stf1(k,l),k=1,3),l=1,3)
c            end if
            if(i.eq.2.or.j.eq.2) then
               a(i,j,1)=0.0d0
            else
               a(i,j,1)=stf1(1,2)+stf2(1,3)
            end if
            if(j.eq.2) then
               a(i,j,2)=0.0d0
            else 
               a(i,j,2)=stf2(2,3)+stf3(1,3)
            end if
            if(i.eq.2) then
               a(i,j,4)=0.0d0
            else
               a(i,j,4)=stf1(2,3)+stf4(1,2)
            end if
            a(i,j,5)=stf1(2,2)+stf2(3,3)+stf3(3,3)
     +              +stf4(2,2)+stf5(1,1)+stf6(1,1)
            if(i.eq.nx) then
               a(i,j,6)=0.0d0
            else
               a(i,j,6)=stf3(2,3)+stf6(1,2)
            end if
            if(j.eq.ny) then
               a(i,j,8)=0.0d0
            else
               a(i,j,8)=stf4(2,3)+stf5(1,3)
            end if
            if(i.eq.nx.or.j.eq.ny) then
               a(i,j,9)=0.0d0
            else
               a(i,j,9)=stf5(1,2)+stf6(1,3)
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
         rhs(1,j)=exact(xleft+0.0d0,yleft+(j-1)*hy)
         rhs(nx+1,j)=exact(xleft+nx*hx,yleft+(j-1)*hy)
      end do

      ! Boundary nodes: top and bottom
      do i=1,nx+1
         rhs(i,1)=exact(xleft+(i-1)*hx,yleft+0.0d0)
         rhs(i,ny+1)=exact(xleft+(i-1)*hx,yleft+ny*hy)
      end do

      ! Interior nodes
      do j=2,ny
         do i=2,nx
            call get_force(value1,2,2,i-1,j-1)
            call get_force(value2,3,1,i-1,j-1)
            call get_force(value3,3,2,i,j-1)
            call get_force(value4,2,1,i-1,j)
            call get_force(value5,1,2,i,j)
            call get_force(value6,1,1,i,j)
            rhs(i,j)=value1+value2+value3
     +              +value4+value5+value6
            if(i.eq.2.and.j.ne.2) then
		value=perm(i-1,j-1)
               call get_stiff(stf1,2,i-1,j-1,value)
		value=perm(i-1,j-1)
               call get_stiff(stf2,1,i-1,j-1,value)
		value=perm(i-1,j)
               call get_stiff(stf4,1,i-1,j,value)
               rhs(i,j)=rhs(i,j)
     +                 -(stf1(1,2)+stf2(1,3))*rhs(i-1,j-1)
     +                 -(stf1(2,3)+stf4(1,2))*rhs(i-1,j)
            end if
            if(i.eq.nx.and.j.ne.ny) then
		value=perm(i,j-1)
               call get_stiff(stf3,2,i,j-1,value)
		value=perm(i,j)
               call get_stiff(stf5,2,i,j,value)
		value=perm(i,j)
               call get_stiff(stf6,1,i,j,value)
               rhs(i,j)=rhs(i,j)
     +                 -(stf3(2,3)+stf6(1,2))*rhs(i+1,j)
     +                 -(stf5(1,2)+stf6(1,3))*rhs(i+1,j+1)
            end if
            if(j.eq.2.and.i.ne.2) then
		value=perm(i-1,j-1)
               call get_stiff(stf1,2,i-1,j-1,value)
		value=perm(i-1,j-1)
               call get_stiff(stf2,1,i-1,j-1,value)
		value=perm(i,j-1)
               call get_stiff(stf3,2,i,j-1,value)
               rhs(i,j)=rhs(i,j)
     +                 -(stf1(1,2)+stf2(1,3))*rhs(i-1,j-1)
     +                 -(stf2(2,3)+stf3(1,3))*rhs(i,j-1)
            end if
            if(j.eq.ny.and.i.ne.nx) then
		value=perm(i-1,j)
               call get_stiff(stf4,1,i-1,j,value)
		value=perm(i,j)
               call get_stiff(stf5,2,i,j,value)
		value=perm(i,j)
               call get_stiff(stf6,1,i,j,value)
               rhs(i,j)=rhs(i,j)
     +                 -(stf4(2,3)+stf5(1,3))*rhs(i,j+1)
     +                 -(stf5(1,2)+stf6(1,3))*rhs(i+1,j+1)
            end if
            if(i.eq.2.and.j.eq.2) then
		value=perm(i-1,j-1)
               call get_stiff(stf1,2,i-1,j-1,value)
		value=perm(i-1,j-1)
               call get_stiff(stf2,1,i-1,j-1,value)
		value=perm(i,j-1)
               call get_stiff(stf3,2,i,j-1,value)
		value=perm(i-1,j)
               call get_stiff(stf4,1,i-1,j,value)
               rhs(i,j)=rhs(i,j)
     +                 -(stf1(1,2)+stf2(1,3))*rhs(i-1,j-1)
     +                 -(stf1(2,3)+stf4(1,2))*rhs(i-1,j)
     +                 -(stf2(2,3)+stf3(1,3))*rhs(i,j-1)
            end if
            if(i.eq.nx.and.j.eq.ny) then
		value=perm(i,j-1)
               call get_stiff(stf3,2,i,j-1,value)
		value=perm(i-1,j)
               call get_stiff(stf4,1,i-1,j,value)
		value=perm(i,j)
               call get_stiff(stf5,2,i,j,value)
		value=perm(i,j)
               call get_stiff(stf6,1,i,j,value)
               rhs(i,j)=rhs(i,j)
     +                 -(stf4(2,3)+stf5(1,3))*rhs(i,j+1)
     +                 -(stf5(1,2)+stf6(1,3))*rhs(i+1,j+1)
     +                 -(stf3(2,3)+stf6(1,2))*rhs(i+1,j)
            end if
         end do
      end do

      ! change to the multigrid solver format
      do j=1,nyf  
         do i=1,nxf
            call assign_mrhs(mrhs,i,j,rhs(i,j),nxf,nyf)
         end do
      end do

      print*, 'go to mgd9v'
c	if (nw.gt.1) iprep=1
      call mgd9v(levels,nxc,nyc,nxf,nyf,nm,iout,istart,iprep,
     +     maxit,tol,mrhs,ma,p,vb,work,ldu,wa,wb,resno,ifail)
      print*,' ifail: ',ifail,' resno: ',resno
      print*, 'out of mgd9v'



c	output the 'exact' solution
c

	open(2,file='../dat/result/exact_2048.dat',status='unknown')
C	write(2,*) nxf,nyf, (p(i),i=1,nxf*nyf)
	write(2,*) nxf-1,nyf-1,nm
           do j=1,nyf
            do i=1,nxf
             write(2,*)p((j-1)*nyf+i)
            enddo
           enddo
	close(2)
        open(1,file='../dat/result/shape2048info.dat',status='unknown')
        write(1,*)nx,ny
        do j=1,ny+1
            write(1,*) 2049
        enddo
        close(1)
c	----------------------------------

      stop
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
      subroutine get_stiff(stf,index,i,j,perm)
      common /param/nx,ny,hx,hy,xleft,yleft,xright,yright
      real*8 hx,hy,stf(3,3),area,x,y,bk1,bk2,bl1,bl2,perm
      real*8 basis_x,basis_y,coor,d,xleft,yleft,xright,yright
      real*8 permc,coef
      integer nx,ny,index,i,j,node,k,l

      area=0.5d0*hx*hy
       permc=0.d0
C.....Compute the barycentric coordinate
c      x=0.0d0
c      y=0.0d0
c      do node=1,3
c         x=x+coor(index,node,1,i,j)
c         y=y+coor(index,node,2,i,j)
c      end do
c      x=x/3.0d0
c      y=y/3.0d0
c      permc=coef(x,y)
      
       do node=1,3
        x=coor(index,node,1,i,j)
        y=coor(index,node,2,i,j)
        permc=permc+coef(x,y)
       enddo
       permc=permc/3.d0

      do k=1,3
         do l=1,3
               bk1=basis_x(index,k,i,j)
               bk2=basis_y(index,k,i,j)
               bl1=basis_x(index,l,i,j)
               bl2=basis_y(index,l,i,j)
               stf(k,l)=area*
     +              permc*(bk1*bl1+bk2*bl2)
         end do
      end do

      return
      end

C     -------------------------------------------------
C     SUBROUTINE GET_FORCE
C     --------------------------------------------------
      ! The integral of the force*basis(node) 
      ! over the element (i,j,index)
      subroutine get_force(value,node,index,i,j)
      common /param/nx,ny,hx,hy,xleft,yleft,xright,yright
      real*8 hx,hy,value,x,y,area,xleft,yleft,xright,yright 
      real*8 force,coor,basis
      integer nx,ny,index,i,j,node

      area=0.5d0*hx*hy

      value=0.0d0
      do k=1,3
         x=coor(index,k,1,i,j)
         y=coor(index,k,2,i,j)
         value=value+force(x,y)*basis(index,node,x,y,i,j)/3.0d0
      end do
      value=value*area

      return
      end

C     -------------------------------------------------
C     SUBROUTINE GET_NUM
C     --------------------------------------------------
      subroutine get_num(num,node,index,i,j)
      common /param/nx,ny,hx,hy,xleft,yleft,xright,yright
      real*8 hx,hy,xleft,yleft,xright,yright 
      integer num,node,index,i,j,nx,ny

      if (index.eq.2) go to 100

      if     (node.eq.1) then
        num=(j-1)*(nx+1)+i
      elseif (node.eq.2) then
        num=(j-1)*(nx+1)+i+1
      elseif (node.eq.3) then
        num=j*(nx+1)+i+1
      else
        print*, 'problem here in get_num!'
      endif
      go to 200

  100 continue

      if     (node.eq.1) then
        num=(j-1)*(nx+1)+i
      elseif (node.eq.2) then
        num=j*(nx+1)+i+1
      elseif (node.eq.3) then
        num=j*(nx+1)+i
      else
        print*, 'problem here in get_num!'
      endif

  200 continue
      return
      end

C     --------------------------------------------
C     FUNCTION COOR
C     --------------------------------------------
      real*8 function coor(index,node,kcoor,i,j)
      common /param/ nx,ny,hx,hy,xleft,yleft,xright,yright
      real*8 hx,hy,xleft,yleft,xright,yright,value1,valu2
      integer nx,ny,index,node,kcoor,i,j

C     print*, 'this is coor'

      if (index.eq.2) go to 100
      
      if     (node.eq.1.and.kcoor.eq.1) then
        coor=xleft+i*hx
      elseif (node.eq.1.and.kcoor.eq.2) then
        coor=yleft+(j-1)*hy+0.5d0*hy
      elseif (node.eq.2.and.kcoor.eq.1) then
        coor=xleft+(i-1)*hx+0.5d0*hx
      elseif (node.eq.2.and.kcoor.eq.2) then
        coor=yleft+(j-1)*hy+0.5d0*hy
      elseif (node.eq.3.and.kcoor.eq.1) then
        coor=xleft+(i-1)*hx+0.5d0*hx
      elseif (node.eq.3.and.kcoor.eq.2) then
        coor=yleft+(j-1)*hy
      else
        print*, 'problem here in function coor!'
      endif
      go to 200
  100 continue 
      
      if     (node.eq.1.and.kcoor.eq.1) then
        coor=xleft+(i-1)*hx+0.5d0*hx
      elseif (node.eq.1.and.kcoor.eq.2) then
        coor=yleft+j*hy
      elseif (node.eq.2.and.kcoor.eq.1) then
        coor=xleft+(i-1)*hx
      elseif (node.eq.2.and.kcoor.eq.2) then
        coor=yleft+(j-1)*hy+0.5d0*hy
      elseif (node.eq.3.and.kcoor.eq.1) then
        coor=xleft+(i-1)*hx+0.5d0*hx
      elseif (node.eq.3.and.kcoor.eq.2) then
        coor=yleft+(j-1)*hy+0.5d0*hy
      else
        print*, 'problem here in function coor!'
      endif

  200 continue

      return
      end

C     ----------------------------------------------
C     FUNCTION BASIS
C     ----------------------------------------------
      real*8 function basis(index,node,x,y,i,j)
      common /param/ nx,ny,hx,hy,xleft,yleft,xright,yright
      real*8 hx,hy,x,y,x0,y0,xleft,yleft,xright,yright 
      integer nx,ny,index,node,kcoor,i,j

      x0=xleft+(i-1)*hx
      y0=yleft+(j-1)*hy

      if (index.eq.2) go to 100

      if     (node.eq.1) then
        basis=-(x-x0-hx)/hx
      elseif (node.eq.2) then
        basis=(x-x0)/hx-(y-y0)/hy
      elseif (node.eq.3) then
        basis=(y-y0)/hy
      else
        print*, 'problem here in function basis!'
      endif
      go to 200

  100 continue

      if     (node.eq.1) then
        basis=-(y-y0-hy)/hy
      elseif (node.eq.2) then
        basis=(x-x0)/hx
      elseif (node.eq.3) then
        basis=(y-y0)/hy-(x-x0)/hx
      else
        print*, 'problem here in function basis!'
      endif


  200 continue

      return
      end

C     ----------------------------------------------
C     FUNCTION BASIS_X
C     ----------------------------------------------
      real*8 function basis_x(index,node,i,j)
      common /param/ nx,ny,hx,hy,xleft,yleft,xright,yright
      real*8 hx,hy,x,y,xleft,yleft,xright,yright 
      integer nx,ny,index,node,kcoor,i,j

      if (index.eq.2) go to 100

      if     (node.eq.1) then
        basis_x=-1.0d0/hx
      elseif (node.eq.2) then
        basis_x=1.0d0/hx
      elseif (node.eq.3) then
        basis_x=0.0d0
      else
        print*, 'problem here in function basis!'
      endif
      go to 200

  100 continue

      if     (node.eq.1) then
        basis_x=0.0d0
      elseif (node.eq.2) then
        basis_x=1.0d0/hx
      elseif (node.eq.3) then
        basis_x=-1.0d0/hx
      else
        print*, 'problem here in function basis!'
      endif

  200 continue

      return
      end

C     ----------------------------------------------
C     FUNCTION BASIS_Y
C     ----------------------------------------------
      real*8 function basis_y(index,node,i,j)
      common /param/ nx,ny,hx,hy,xleft,yleft,xright,yright
      real*8 hx,hy,x,y,xleft,yleft,xright,yright 
      integer nx,ny,index,node,kcoor,i,j

      if (index.eq.2) go to 100

      if     (node.eq.1) then
        basis_y=-0.0d0
      elseif (node.eq.2) then
        basis_y=-1.0d0/hy
      elseif (node.eq.3) then
        basis_y=1.0d0/hy
      else
        print*, 'problem here in function basis!'
      endif
      go to 200

  100 continue

      if     (node.eq.1) then
        basis_y=-1.0d0/hy
      elseif (node.eq.2) then
        basis_y=0.0d0
      elseif (node.eq.3) then
        basis_y=1.0d0/hy
      else
        print*, 'problem here in function basis!'
      endif

  200 continue

      return
      end


C*****************************************************
C        INPUT DATA
C*****************************************************

C     -------------------------------------------
C     SOURCE FUNCTION
C     -------------------------------------------
      real*8 function force(x,y)
      real*8 x,y,pi,eps
      common/epsilon/eps
       pi=dacos(-1.0d0)
      force=-1.0d0   

      return
      end

C     -------------------------------------------
C     EXACT SOLUTION
C     -------------------------------------------
      real*8 function exact(x,y)
      real*8 x,y,pi,eps
       common/epsilon/eps
       pi=dacos(-1.0d0)
       exact=0.0d0
       return
      end

