C       MULTISCALE FINITE ELEMENT METHOD           *
C       OVERSAMPLING                               *
C       BOUNDARY CONDITION: Dirichlet              *
C       TO SOLVE THE MATRIX PROBLEM, USE dmgd9v    *
C                                                  *
C***************************************************
      program main
      parameter (max_nm=1402193, ! nm for level 9
     +           max_nxf=1025,  ! nxf for level 9
     +           max_nyf=1025)  ! nyf for level 9

      real*8 ma(max_nm*9),mrhs(max_nm),a(max_nxf,max_nyf,9),
     +       vb(max_nm),ldu(3*max_nm),work(max_nxf*12),
     +       wa(max_nm),wb(max_nm),p(max_nm),rhs(max_nxf,max_nyf),
     +	     perm(-511:2560,-511:2560),
     +        b(2,3,max_nm)
      real*8 stiff(4,4,1024,1024) ! the size of global meshes
                             ! nxglobal x nyglobal
C	real*8 rhs1(3,3),rhs2(3,3)
      real*8 rhs_elem(4,1024,1024),rhs1(4)
      real*8 hx,hy,tol,resno,value,area,x,y,eu,eflux
      real*8 stf1(4,4),stf2(4,4),stf3(4,4),stf4(4,4),stf(4,4)
      real*8 value1,value2,value3,value4,eps
      real*8 coor,basis,exact,flux,basis_x,basis_y
      real*8 hxl,hyl,xl,yl,b1,b2,b3,local_basis_x,local_basis_y
      integer nx,ny,nxc,nyc,nxf,nyf,nm,levels,i,j,k,l,
     +        istart,iprep,maxit,ifail,index,kcoor,node,
     +        iout(6),n1,n2,n3,num,mil,mjl,nb,ntran,mmi,mmj
	real*8 xleft,yleft,xright,yright
	integer nxfine,nyfine
	real*8 p_ex(5600594),p_ex_c(1050625),p_max,norml2,v1
         character*13 filename

        common /epsilon/eps

      common /param/ nx,ny,hx,hy,xleft,yleft

      data nxc,nyc,istart,iprep,maxit,tol,ifail
     +    /5,5,0,0,50,1.0d-8,110/
      data iout/0,0,0,0,0,0/

C.....PARAMETERS..............................

      open(1, file='input.dat',status='old')
      read(1,*) levelslocal,nxlocal,nxglobal,levelsglobal,ineps,
     +          mimax,mjmax
	read(1,*)
	read(1,*) xleft,yleft,xright,yright
c	read(1,*)
	eps=1.d0/ineps

      close(1)
        pi=dacos(-1.d0)

      levels=levelsglobal     ! multigrid levels
      nx=nxglobal        ! number of spaces in x direction
      ny=nx        ! number of spaces in y direction

      hx=(xright-xleft)/nx
      hy=(yright-yleft)/ny
      area=hx*hy   ! area of the element

C.....LOCAL PARAMETERS.........................

      nxl=nxlocal
      nyl=nxl
      hxl=hx/nxl
      hyl=hy/nyl

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


c*********************************************
c       input the permeability data
c*********************************************
c       nxstart=-255
c       nystart=-255
c       open(1,file='../dat/perms.dat')
c           read(1,*) nx1,ny1,ntot1
c           do i=nxstart,nx1+nxstart-1
c               read(1,*) (perm(i,j),j=nystart,ny1+nystart-1)
c           enddo
c       close(1)

c     do i=1,nx*nxl
c        do j=1,ny*nyl
c           perm(i,j)=(coef((i-1)*hxl,(j-1)*hyl)+coef((i-1)*hxl,j*hyl)
c    +                +coef(i*hxl,(j-1)*hyl)+coef(i*hxl,j*hyl))/4.0d0
c        end do
c     end do
   


C.....READ THE ELEMENT STIFFNESS MATRIX.........

      open(1,file='base/stf',status='old')
      open(2,file='base/rhs',status='old')
      do 100 mmj=1,mjmax
      do 100 mmi=1,mimax
         read(1,*) ((stf(i,j),j=1,4),i=1,4)   
         read(2,*) (rhs1(i),i=1,4) 
         do mj=mmj,ny,mjmax
         do mi=mmi,nx,mimax
            do k=1,4
            do l=1,4
               stiff(k,l,mi,mj)=stf(k,l)
            end do
               rhs_elem(k,mi,mj)=rhs1(k)
            end do
         end do
         end do
100	continue

      close(1)
      close(2)
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
            call get_stiff(stf1,i-1,j-1,stiff)
            call get_stiff(stf2,i,j-1,stiff)
            call get_stiff(stf3,i-1,j,stiff)
            call get_stiff(stf4,i,j,stiff)
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
         rhs(1,j)=exact(xleft,yleft+(j-1)*hy)
         rhs(nx+1,j)=exact(xleft+nx*hx,yleft+(j-1)*hy)
      end do

      ! Boundary nodes: top and bottom
      do i=1,nx+1
         rhs(i,1)=exact(xleft+(i-1)*hx,yleft)
         rhs(i,ny+1)=exact(xleft+(i-1)*hx,yleft+ny*hy)
      end do 

      ! Interior nodes
      do j=2,ny
         do i=2,nx
c            call get_force(value1,3,i-1,j-1,rhs_elem)
c            call get_force(value2,4,i,j-1,rhs_elem)
c            call get_force(value3,2,i-1,j,rhs_elem)
c            call get_force(value4,1,i,j,rhs_elem)
            value1=rhs_elem(3,i-1,j-1)
            value2=rhs_elem(4,i,j-1)
            value3=rhs_elem(2,i-1,j)
            value4=rhs_elem(1,i,j)
            rhs(i,j)=value1+value2+value3+value4
            if(i.eq.2.and.j.ne.2.and.j.ne.ny) then
               call get_stiff(stf1,i-1,j-1,stiff)
               call get_stiff(stf3,i-1,j,stiff)
               rhs(i,j)=rhs(i,j)
     +            -stf1(1,3)*rhs(i-1,j-1)
     +            -(stf1(4,3)+stf3(1,2))*rhs(i-1,j)
     +            -stf3(4,2)*rhs(i-1,j+1)
            end if
            if(i.eq.nx.and.j.ne.2.and.j.ne.ny) then
               call get_stiff(stf2,i,j-1,stiff)
               call get_stiff(stf4,i,j,stiff)
               rhs(i,j)=rhs(i,j)
     +                 -stf2(2,4)*rhs(i+1,j-1)
     +                 -(stf2(3,4)+stf4(2,1))*rhs(i+1,j)
     +                 -stf4(3,1)*rhs(i+1,j+1)
            end if
            if(j.eq.2.and.i.ne.2.and.i.ne.nx) then
               call get_stiff(stf1,i-1,j-1,stiff)
               call get_stiff(stf2,i,j-1,stiff)
               rhs(i,j)=rhs(i,j)
     +            -stf1(1,3)*rhs(i-1,j-1)
     +            -(stf1(2,3)+stf2(1,4))*rhs(i,j-1)
     +            -stf2(2,4)*rhs(i+1,j-1)
            end if
            if(j.eq.ny.and.i.ne.2.and.i.ne.nx) then
               call get_stiff(stf3,i-1,j,stiff)
               call get_stiff(stf4,i,j,stiff)
               rhs(i,j)=rhs(i,j)
     +            -stf3(4,2)*rhs(i-1,j+1)
     +            -(stf3(3,2)+stf4(4,1))*rhs(i,j+1)
     +            -stf4(3,1)*rhs(i+1,j+1)
            end if
            if(i.eq.2.and.j.eq.2) then
               call get_stiff(stf1,i-1,j-1,stiff)
               call get_stiff(stf2,i,j-1,stiff)
               call get_stiff(stf3,i-1,j,stiff)
               rhs(i,j)=rhs(i,j)
     +            -stf1(1,3)*rhs(i-1,j-1)
     +            -(stf1(2,3)+stf2(1,4))*rhs(i,j-1)
     +            -stf2(2,4)*rhs(i+1,j-1)
     +            -(stf1(4,3)+stf3(1,2))*rhs(i-1,j)
     +            -stf3(4,2)*rhs(i-1,j+1)
            end if
            if(i.eq.2.and.j.eq.ny) then
               call get_stiff(stf1,i-1,j-1,stiff)
               call get_stiff(stf3,i-1,j,stiff)
               call get_stiff(stf4,i,j,stiff)
               rhs(i,j)=rhs(i,j)
     +            -stf1(1,3)*rhs(i-1,j-1)
     +            -(stf1(4,3)+stf3(1,2))*rhs(i-1,j)
     +            -stf3(4,2)*rhs(i-1,j+1)
     +            -(stf3(3,2)+stf4(4,1))*rhs(i,j+1)
     +            -stf4(3,1)*rhs(i+1,j+1)
            end if              
            if(i.eq.nx.and.j.eq.2) then
               call get_stiff(stf1,i-1,j-1,stiff)
               call get_stiff(stf2,i,j-1,stiff)
               call get_stiff(stf4,i,j,stiff)
               rhs(i,j)=rhs(i,j)
     +            -stf1(1,3)*rhs(i-1,j-1)
     +            -(stf1(2,3)+stf2(1,4))*rhs(i,j-1)
     +            -stf2(2,4)*rhs(i+1,j-1)
     +            -(stf2(3,4)+stf4(2,1))*rhs(i+1,j)
     +            -stf4(3,1)*rhs(i+1,j+1)
            end if
            if(i.eq.nx.and.j.eq.ny) then
               call get_stiff(stf2,i,j-1,stiff)
               call get_stiff(stf3,i-1,j,stiff)
               call get_stiff(stf4,i,j,stiff)
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

      print*, 'go to mgd9v'
      call mgd9v(levels,nxc,nyc,nxf,nyf,nm,iout,istart,iprep,
     +     maxit,tol,mrhs,ma,p,vb,work,ldu,wa,wb,resno,ifail)
      print*,' ifail: ',ifail,' resno: ',resno
      print*, 'out of mgd9v'

c     print*, 'p=',(p(i),i=1,nxf*2)
c	output the result
	open(2,file='result/Omsfem_solution.dat')
	write(2,*) nxf,nyf,(p(i),i=1,nxf*nyf)
	close(2)

c	compare with the exact solution

 	open(1,file='../dat/result/exact_2048.dat',
     +          form='formatted')
  	read(1,*) nxfine,nyfine,(p_ex(i),i=1,nxfine*nyfine)
 	close(1)

c.	average on the fine grid to yiled coarse grid 'exact' solution

c       numx=(nxfine-1)/nx 
c       numy=(nyfine-1)/ny 
c       do j=1,ny 
c       do i=1,nx 
c          v1=0.d0
c          k0=(j-1)*numy*nxfine+(i-1)*numx+1
c          do jj=1,numy+1
c          do ii=1,numx+1
c               kk=k0+(jj-1)*nxfine+ii-1
c               v1=v1+p_ex(kk)
c          enddo
c          enddo
c          v1=v1/(numx+1)/(numy+1)
c          k=(j-1)*nx+i
c          p_ex_c(k)=v1
c       enddo
c       enddo

	eu=0.d0
	value=0.d0
        norml2=0.d0
        p_max=0.0d0
	nxfine=nxfine-1
	nyfine=nyfine-1
 	do j=2,ny
 	   do i=2,nx
 		n1=(j-1)*(nx+1)+i
 		n2=(j-1)*(nyfine/ny)*(nxfine+1)+(i-1)*(nxfine/nx)+1
 		eu=eu+(p(n1)-p_ex(n2))**2
 		if(value.lt.dabs(p(n1)-p_ex(n2))) value=dabs(p(n1)-p_ex(n2))
                norml2=norml2+p_ex(n2)**2
                if(p_max.lt.dabs(p_ex(n2))) p_max=dabs(p_ex(n2))
c		eu=eu+(p(n1)-p_ex_c(n1))**2
c		if(value.lt.abs(p(n1)-p_ex_c(n1))) value=abs(p(n1)-p_ex_c(n1))
c               norml2=norml2+p_ex_c(n1)**2
c               if(p_max.lt.abs(p_ex(n1))) p_max=abs(p_ex_c(n1))
 	   enddo
 	enddo
        norml2=dsqrt(norml2*hx*hy)  
 	eu=dsqrt(eu*hx*hy)
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

            filename = 'error'//char(m1+48)//char(48+m2)//char(48+m3)//
     $           char(48+m4)//'.dat'

 
        open(12,file=filename)
 	write(12,*) " nx,  relative L^2-error,  relative  L^\\infty-error"
 	write(12,*) nx,eu/norml2,value/p_max
 	close(12)
 	print*, 'L^2 error is  ', eu,eu/norml2
 	print*, "L^\\infty  error is  ",  value,value/p_max
        print*, 'l2-norm and L_inf norm of exact solution: ', norml2,
     &  p_max
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
C----------------------------------------------------
C     SUBROUTINE GET_STIFF
C----------------------------------------------------
      subroutine get_stiff(stf,i,j,stiff)
      real*8 stf(4,4),stiff(4,4,1024,1024)
      integer i,j,k,l,nx,ny

      do 20 k=1,4
         do 10 l=1,4
         stf(k,l)=stiff(k,l,i,j)
 10      continue
 20   continue
 
      return
      end

C     -------------------------------------------------
C     SUBROUTINE GET_FORCE
C     --------------------------------------------------
      ! The integral of the force*basis(node) 
      ! over the element (i,j)
      subroutine get_force(value,node,i,j,rhs_elem)
      common /param/nx,ny,hx,hy,xleft,yleft
      real*8 hx,hy,value,x,y,area,xleft,yleft
      real*8 force,coor_gauss,basis,rhs_elem(4,1024,1024)
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
C     --------------------------------------------
C     FUNCTION COOR
C     --------------------------------------------
      real*8 function coor(node,kcoor,i,j)
      common /param/ nx,ny,hx,hy,xleft,yleft
      real*8 hx,hy,x0,y0,xleft,yleft
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
      common /param/ nx,ny,hx,hy,xleft,yleft
      real*8 hx,hy,x0,y0,tmp,xleft,yleft
      integer nx,ny,index,node,kcoor,i,j

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
      common /param/ nx,ny,hx,hy,xleft,yleft
      real*8 hx,hy,x,y,x0,y0,area,xleft,yleft
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

C     -------------------------------------------
C     EXACT SOLUTION
C     -------------------------------------------
      real*8 function exact(x,y)
      real*8 x,y
      exact=0.0d0
c     exact=1.0d0-x
      end
