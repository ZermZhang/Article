C 2012年03月13日 星期二 17时56分29秒 2012年03月17日 星期六 13时51分43秒 
C  *************************************************
C     Compare at the finest mesh using the rule three 
C     midpoints of the edges
C***************************************************
      
      subroutine get_error(lgmap,nel,mxmax,i_fnum,mxy,iflag,
     &            u,kdof,np,nz,p,symm,sigma1,sigma2,ineps)
      implicit none
      integer nel,mxy,kdof,np,nz,ineps
      double precision symm,sigma1,sigma2
      integer i_fnum(mxy),mxmax(my),lgmap(nel,4),
     &   iflag(-1:mx+1,-1:my+1)
      double precision p(4,np,np,nz),u(kdof)
      
      double precision xmax,ymax,hx,hy,xleft,yleft,eps
      integer nx,ny,mx,my
      common /domain/mx,my,xmax,ymax,xleft,yleft  !粗网格尺寸，参数
      common /param/nx,ny,hx,hy !细网格尺寸，参数
      common/epsilon/eps

      double precision uhf(4),uhm(4),uec(4),bk1(4),bk2(4)     
      double precision,allocatable :: u_ex(:)
      integer,allocatable :: mxmax1(:)
      character*35 filename      
      double precision eu,value,enorml2,p_max,h1error,h1exact,pm,utmp,
     +                   uhtmp,dctmp,basis_x,basis_y,xf,yf
      integer nxfine,nyfine,nm2,k,mi,mj,nxl,nyl,mindex,i,j,ni,nj,k1,k2,
     +         nindex,indexl,node,nodel,il,jl,nl,iter,ii1,ii2,ii3,m1,m2,
     +         m3,m4,m,ilf,jlf
      double precision hxf,hyf,area,el2tmp,eutmp,h1tmp,euxtmp,basisf
 
      eu=0.d0
      value=0.d0
      enorml2=0.d0
      p_max=0.0d0
      h1error=0.d0
      h1exact=0.d0
C-------------------------------------------------
C     compare with the exact solution
      open(1,file='../dat/shape2048info.dat')
c      open(1,file='../dat/shape1024info.dat')
      read(1,*)nxfine,nyfine
      allocate(mxmax1(1+nyfine))
      do j=1,nyfine+1
       read(1,*)mxmax1(j)
      enddo
      close(1)

      k=0
c      open(2,file='../dat/linear_1024.dat')
      open(2,file='../dat/exact_2048.dat')
      read(2,*)nxfine,nyfine,nm2
      allocate(u_ex(nm2))
      do 551 mj=1,nyfine+1
      do 551 mi=1,mxmax1(mj)
        k=k+1
        read(2,*)u_ex(k)
551   continue
      close(2)
      print*,k,nm2,u_ex(130303)
      
      nxl=nxfine/mx/nx
      nyl=nyfine/my/ny
      print*,nxfine,mx,nx,nxl
      
      hxf=hx/nxl
      hyf=hy/nyl
      area=hxf*hyf   ! area of the finest grid
C-------------------------------------------------
      do 3000 mj=1,my
      do 3000 mi=1,mxmax(mj)
       if(iflag(mi,mj).eq.0)then
         call getel_num(k2,mxmax,i_fnum,mxy,mi,mj,0,0,0)
         do 180 i=1,nx  
          ni=(mi-1)*nx+i    
          do 170 j=1,ny
c           print*,'i,j',i,j
           nj=(mj-1)*ny+j
            do node=1,4
	      call get_num(k1,node,i,j)
             utmp=0.d0
             do m=1,4
              pm=p(m,mi,mj,k1)!basis(mindex,m,x,y,mi,mj)
              utmp=utmp+u(lgmap(k2,m))*pm
             enddo
             uhm(node)=utmp
            enddo
            do 1801 il=1,nxl    
            do 1701 jl=1,nyl
C             print*,'il,jl,indexl',il,jl,indexl
            do nodel=1,4
             call get_numf(ilf,jlf,nodel,mi,mj,i,j,
     &              il,jl,nxl,nyl)
             xf=xleft+(ilf-1)*hxf
             yf=yleft+(jlf-1)*hyf
             uhtmp=0.d0
             do m=1,4
              uhtmp=uhtmp+uhm(m)*basisf(m,xf,yf,ni,nj)
             enddo
             uhf(nodel)=uhtmp
             call getnumf(nl,mxmax1,ilf,jlf,nxfine,nyfine)
             uec(nodel)=u_ex(nl)        
             dctmp=dabs(uhf(nodel)-uec(nodel))
             if(value.lt.dctmp) value=dctmp
             if(p_max.lt.dabs(uec(nodel))) p_max=dabs(uec(nodel))
C             bk1(nodel)=basis_x(nodel,hxf,hyf,xf,yf)
C             bk2(nodel)=basis_y(nodel,hxf,hyf,xf,yf)
            enddo
            call get_integral(el2tmp,eutmp,h1tmp,euxtmp,uhf,uec,
     &          hxf,hyf)
            eu=eu+el2tmp
            enorml2=enorml2+eutmp
            h1error=h1error+h1tmp
            h1exact=h1exact+euxtmp
1701       continue
1801      continue
160      continue
170      continue
180     continue
       else
        do j=1,ny
        do i=1,nx
         ni=(mi-1)*nx+i 
         nj=(mj-1)*ny+j
         nindex=1 
          call getel_num(k2,mxmax,i_fnum,mxy,mi,mj,nindex,i,j)
           do 80 il=1,nxl    
           do 70 jl=1,nyl
            do nodel=1,4
             call get_numf(ilf,jlf,nodel,mi,mj,i,j,
     &              il,jl,nxl,nyl)
             xf=xleft+(ilf-1)*hxf
             yf=yleft+(jlf-1)*hyf
             uhtmp=0.d0
             do m=1,4
              uhtmp=uhtmp+u(lgmap(k2,m))*basisf(m,xf,yf,ni,nj)
             enddo
             uhf(nodel)=uhtmp
             call getnumf(nl,mxmax1,ilf,jlf,nxfine,nyfine)
             uec(nodel)=u_ex(nl)        
             dctmp=dabs(uhf(nodel)-uec(nodel))
             if(value.lt.dctmp) value=dctmp
             if(p_max.lt.dabs(uec(nodel))) p_max=dabs(uec(nodel))
C              bk1(nodel)=basis_x(indexl,nodel,hxf,hyf)
C              bk2(nodel)=basis_y(indexl,nodel,hxf,hyf)
            enddo
            call get_integral(el2tmp,eutmp,h1tmp,euxtmp,uhf,uec,
     &          hxf,hyf)
            eu=eu+el2tmp
            enorml2=enorml2+eutmp
            h1error=h1error+h1tmp
            h1exact=h1exact+euxtmp
70       continue
80      continue

        enddo
        enddo
       endif

 3000     continue

C-------------------------------------------------
C        enorml2=dsqrt(enorml2*area)
C  	 eu=dsqrt(eu*area)
C        h1error=dsqrt(h1error*area)
C        h1exact=dsqrt(h1exact*area)
	 enorml2=dsqrt(enorml2)
 	 eu=dsqrt(eu)
       h1error=dsqrt(h1error)
       h1exact=dsqrt(h1exact)
cccccccccccccccccccccccccc
c  This creates a file called error????.dat, where ???? is the partition
c  number on x-direction.
cccccccccccccccccccccccccc
       iter=int(1.0/eps)
            ii1 = iter / 1000
            m1 = mod(ii1,10)
            ii2 = iter / 100
            m2 = mod(ii2,10)
            ii3 = iter / 10
            m3 = mod(ii3,10)
            m4 = mod(iter,10)

       filename = '../dat/result/j1nciperror'//char(m1+48)//char(48+m2)
     $          //char(48+m3)//char(48+m4)//'.dat'

        print*,'L^2 error is  ', eu,eu/enorml2

 	 print*, "L^\\infty  error is  ",  value,value/p_max
 	 print*,"H^1  error is  ",  h1error,h1error/h1exact
 	 print*, 'l2, L_inf and H1-norm norm of exact solution: ',
     &       enorml2, p_max,h1exact
        open(12,file=filename)
 	write(12,*) 'eps is  ', eps
 	write(12,*) 'coarse mesh size ', mx
 	write(12,*) 'fine mesh size ', nx
 	write(12,*) 'symm is (-1) SIPDG ', symm
 	write(12,*) 'sigma1=%f,sigma2=f% ', sigma1,sigma2
 	write(12,*) 'L^2 error is  ', eu,eu/enorml2
 	write(12,*) "L^\\infty  error is  ",  value,value/p_max
 	write(12,*) "H^1  error is  ",  h1error,h1error/h1exact
        write(12,*) 'l2, L_inf and H1-norm norm of exact solution: ',
     &       enorml2, p_max,h1exact
 	close(12)
      stop
      end
C     -------------------------------------------------
C     SUBROUTINE integral on the finest grid using the rule of three middle points
C     --------------------------------------------------
      subroutine get_integral_old(el2tmp,eutmp,h1tmp,euxtmp,uhf,uec,
     &          bk1,bk2)
      implicit none
      double precision el2tmp,eutmp,h1tmp,euxtmp
      double precision uhf(4),uec(4),bk1(4),bk2(4)
      double precision h1tmp1,euxtmp1,h1tmp2,euxtmp2
      integer m,l

      el2tmp=0.d0
      eutmp=0.d0
      h1tmp1=0.d0
      euxtmp1=0.d0
      h1tmp2=0.d0
      euxtmp2=0.d0
      do m=1,4
       l=m+1
       if(l.eq.5)l=1
       el2tmp=el2tmp+((uhf(m)-uec(m)+uhf(l)-uec(l))/2.d0)**2
       eutmp=eutmp+((uec(m)+uec(l))/2.d0)**2
       h1tmp1=h1tmp1+(uhf(m)-uec(m))*bk1(m)
       euxtmp1=euxtmp1+uec(m)*bk1(m)
       h1tmp2=h1tmp2+(uhf(m)-uec(m))*bk2(m)
       euxtmp2=euxtmp2+uec(m)*bk2(m)
      enddo
       el2tmp=el2tmp/4.d0            
       eutmp=eutmp/4.d0
       h1tmp=h1tmp1**2+h1tmp2**2
       euxtmp=euxtmp1**2+euxtmp2**2
      end
	
	
	
      subroutine get_integral(el2tmp,eutmp,h1tmp,euxtmp,uhf,uec,
     &          hxf,hyf)
      implicit none
      double precision el2tmp,eutmp,h1tmp,euxtmp,hxf,hyf,area
      double precision uhf(4),uec(4),bk1(4),bk2(4),tmp(4),A(4,4)
      double precision h1tmp1,euxtmp1,h1tmp2,euxtmp2
      integer m,l,i,j
	
	
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
c      tmp=h1tmp*area
      end 
	
	

C     -------------------------------------------------
C     SUBROUTINE GET_NUM 获得精确解细网格层的（ni,nj)编号
C     --------------------------------------------------
      subroutine get_numf(ni,nj,node,mi,mj,i1,j1,
     &                il,jl,nxl,nyl)
      implicit none
      integer ni,nj,node,index,mi,mj,i1,j1,il,jl,nxl,nyl
      
      double precision hx,hy
      integer nx,ny   
      common /param/nx,ny,hx,hy

      integer i,j 
      
      i=((mi-1)*nx+i1-1)*nxl+il
      j=((mj-1)*ny+j1-1)*nyl+jl  !Global fine mesh location


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
C**************************************
C     get the indice of dof for the exact solution 
C
C**************************************
      subroutine getnumf(k,mxmax,mi,mj,nxfine,nyfine)
      implicit none
      integer k,mi,mj,nxfine,nyfine
      integer mxmax(nyfine+1)
      integer j
      k=0
      do j=1,mj-1
        k=k+mxmax(j)
      enddo
      k=k+mi
      return
      end     
C     ----------------------------------------------
C     FUNCTION BASIS 细网格基函数
C     ----------------------------------------------
      function basisf(node,x,y,i,j)
      implicit none
      double precision basisf,x,y
      integer index,node,i,j
      
      double precision xmax,ymax,hx,hy,xleft,yleft,area
      integer nx,ny,mx,my
      common /domain/mx,my,xmax,ymax,xleft,yleft
      common /param/nx,ny,hx,hy

      double precision x1,y1
      
      x1=xleft+(i-1)*hx
      y1=yleft+(j-1)*hy
	area=hx*hy
	
      if     (node.eq.1.) then
        basisf=(x-x1-hx)*(y-y1-hy)/area
      elseif (node.eq.2) then
        basisf=-(x-x1)*(y-y1-hy)/area
      elseif (node.eq.3) then
        basisf=(x-x1)*(y-y1)/area
	elseif (node.eq.4) then
	  basisf=-(x-x1-hx)*(y-y1)/area
      else
        print*, 'problem here in function basis!'
      endif

      return
      end      

