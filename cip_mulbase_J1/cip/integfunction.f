C  2012年03月17日 星期六 13时52分09秒    
C   本程序包包含所有积分程序，程序名含义
C   r_ 右端； p_ 加罚项； element_ 单元积分； edge_ 单元边上积分;
C   f 细网格积分； c 粗网格积分；  fc 悬点情形 细网格对粗网格  cf 粗对细
C****************************************************
C       ELEMENT INFORMATION
C****************************************************
C     FUNCTION Element_Integral
C   单元积分 (perm\nabla \Phi_x^k,\nabla \Phi_x^l)_K
C   积分方法，粗网格上的积分利用粗网格内所有细网格上积分之和
C	have clasified
C     -----------------------------------------------
      subroutine element_integ_c(stfi,mi,mj,np,nz,p)
      implicit none
      integer mindex,mi,mj,np,nz
      double precision p(4,np,np,nz)
      double precision stf(4,4),stfi(4,4),tmp
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft

      integer i,j,index,n1,n2,n3,n4,m,l
      double precision value
      
	
	stfi = 0.0d0
      do 80 i = 1,nx
          do 70 j=1,ny
	    
	        call get_stiff(stf,mi,mj,i,j)
            call get_num(n1,1,i,j)
            call get_num(n2,2,i,j)
            call get_num(n3,3,i,j)
            call get_num(n4,4,i,j)
            do m=1,4
		      do l=1,4
C 			    stfi(m,l)=stfi(m,l)
     			    tmp=p(m,mi,mj,n1)*p(l,mi,mj,n1)*stf(1,1)
     +			       +p(m,mi,mj,n1)*p(l,mi,mj,n2)*stf(1,2)
     +			       +p(m,mi,mj,n1)*p(l,mi,mj,n3)*stf(1,3)
     +			       +p(m,mi,mj,n1)*p(l,mi,mj,n4)*stf(1,4)
     +			       +p(m,mi,mj,n2)*p(l,mi,mj,n1)*stf(2,1)
     +			       +p(m,mi,mj,n2)*p(l,mi,mj,n2)*stf(2,2)
     +			       +p(m,mi,mj,n2)*p(l,mi,mj,n3)*stf(2,3)
     +			       +p(m,mi,mj,n2)*p(l,mi,mj,n4)*stf(2,4)
     +			       +p(m,mi,mj,n3)*p(l,mi,mj,n1)*stf(3,1)
     +			       +p(m,mi,mj,n3)*p(l,mi,mj,n2)*stf(3,2)
     +			       +p(m,mi,mj,n3)*p(l,mi,mj,n3)*stf(3,3)
     +			       +p(m,mi,mj,n3)*p(l,mi,mj,n4)*stf(3,4)
     +			       +p(m,mi,mj,n4)*p(l,mi,mj,n1)*stf(4,1)
     +			       +p(m,mi,mj,n4)*p(l,mi,mj,n2)*stf(4,2)
     +			       +p(m,mi,mj,n4)*p(l,mi,mj,n3)*stf(4,3)
     +			       +p(m,mi,mj,n4)*p(l,mi,mj,n4)*stf(4,4)
			     stfi(m,l)=stfi(m,l)+tmp
            end do
          end do
		  
70        end do
80    end do

      return
      end
C     -----------------------------------------------
C     FUNCTION force right Element_Integral
C   单元积分 (f， \Phi_x^l)_K
C   积分方法，粗网格上的积分利用粗网格内所有细网格上积分之和
C	have clasified
C     -----------------------------------------------
      function r_element_integ_c(mi,mj,m,np,nz,p)
      implicit none
      integer mi,mj,m,np,nz
      double precision p(4,np,np,nz),r_element_integ_c
      double precision stf(4,4)
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft

      double precision volume,f1,f2,f3,f4,force
      integer i,j,ii1,jj1,n1,n2,n3,n4,ntran,index
      
      r_element_integ_c=0.d0
      volume=hx*hx/4.d0
       do 80 i=1,nx
         do 70 j=1,ny
           ii1=ntran((mi-1)*nx+i)
           jj1=ntran((mj-1)*ny+j)
           f1=force(ii1,jj1)
           f2=force(ii1+1,jj1)
           f3=force(ii1+1,jj1+1)
	       f4=force(ii1,jj1+1)
           call get_num(n1,1,i,j)
           call get_num(n2,2,i,j)
           call get_num(n3,3,i,j)
	       call get_num(n4,4,i,j)
           r_element_integ_c=r_element_integ_c
     +           +p(m,mi,mj,n1)*f1*volume
     +           +p(m,mi,mj,n2)*f2*volume
     +           +p(m,mi,mj,n3)*f3*volume
     +           +p(m,mi,mj,n4)*f4*volume
 70      continue
 80    continue


      return
      end
C     -----------------------------------------------
C     FUNCTION Edge_Integral edge parameter 
C     相邻粗网格单元K1，K2在公共单元边e上的积分，(c\nabla\phi_1*n_1, \phi_2)_e
C     \phi_1为单元K1中基函数， \phi_2为单元K2中基函数， n_1 为单位外法向量
C     (nvx，nvy）----边的外法向,如下图表示
C               (0,1)
C         _____________
C         |           |
C         |           |
C         |           |
C   (-1,0)|           |(1,0)
C         |           |
C         |           | 
C         |           |
C         |           |
C         |           |
C         |           |
C         -------------
C            (0,-1)
C     （index1,mi1,mj1） 刻画单元1，（kx1,ky1）记录边在单元1中的位置信息 取 nx-1 or ny-1 or 0
C     （index2,mi2,mj2） 刻画单元1，（kx2,ky2）记录边在单元2中的位置信息 取 nx-1 or ny-1 or 0
C      （li1,lj1)积分边的起点对应单元的位置. 
C	have clasified
C     -----------------------------------------------
      subroutine edge_integ_c(stf,mi1,mj1,kx1,ky1,
     +                mi2,mj2,kx2,ky2,
     +                li1,lj1,nvx,nvy,np,nz,p)
      implicit none
      integer index1,mi1,mj1,kx1,ky1,index2,mi2,mj2,kx2,ky2,
     +         li1,lj1,nvx,nvy,np,nz
      double precision p(4,np,np,nz),stf(4,4)
      
      double precision xmax,ymax,hx,hy,xleft,yleft,x1,x2,y1,y2,
     +                 x0,y0,x3,y3
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft

      double precision xyn2,xn1,yn1,perm1,perm2,value2,value4,
     +           basis_x,basis_y,bx1,by1,bx2,by2,bx3,by3,coef,
     +           value21,value22,value23,value41,value42,value43,
     +           perm3
      integer kx,ky,node2,node4,i,j,i1,j1,i2,j2,ii1,jj1,
     +         ntran,m,l,k,n1,n2,n4,node
     
        stf=0.d0
      
      x0=xleft+dble(mi1-1)*xmax
      y0=yleft+dble(mj1-1)*ymax
	
      xyn2=dsqrt(dble(nvx)**2+dble(nvy)**2)
      xn1=dble(nvx)/xyn2  !单位法向量
      yn1=dble(nvy)/xyn2  !normal direction
      kx=abs(nvy)
      ky=abs(nvx)
      
      value21=hx*(1.0d0  )
      value22=hx*(0.0d0  )
      value23=hx*(1/2.0d0)
      
      value41=hx*(0.0d0  )
      value42=hx*(1.0d0  )
      value43=hx*(1/2.0d0)

      if (xn1 .eq. -1 .and. yn1 .eq. 0) then
	     node2=1
	     node4=4
      else if (xn1 .eq. 1 .and. yn1 .eq. 0) then
	     node2=2
	     node4=3
      else if (xn1 .eq. 0 .and. yn1 .eq. -1) then
         node2=1
         node4=2
      else if  (xn1 .eq. 0 .and. yn1.eq. 1) then
         node2=4
         node4=3
      endif
      do 80 k=1,nx
          i=max(kx*k,1)
          j=max(ky*k,1)
          i1=i+kx1
          j1=j+ky1
          i2=i+kx2
          j2=j+ky2
          ii1=ntran((li1-1)*nx+i1)
          jj1=ntran((lj1-1)*ny+j1)
	      if (xn1 .eq. -1 .and. yn1 .eq. 0) then
	         x1=x0+(i1-1)*hx
	         y1=y0+(j1-1)*hy
		     x2=x0+(i1-1)*hx
		     y2=y0+(j1  )*hy
          else if (xn1 .eq. 1 .and. yn1 .eq. 0) then
	         x1=x0+(i1  )*hx
	         y1=y0+(j1-1)*hy
	         x2=x0+(i1  )*hx
		     y2=y0+(j1  )*hy
          else if (xn1 .eq. 0 .and. yn1 .eq. -1) then
             x1=x0+(i1-1)*hx
	         y1=y0+(j1-1)*hy
             x2=x0+(i1  )*hx
		     y2=y0+(j1-1)*hy
          else if  (xn1 .eq. 0 .and. yn1.eq. 1) then
             x1=x0+(i1-1)*hx
	         y1=y0+(j1  )*hy
             x2=x0+(i1  )*hx
		     y2=y0+(j1  )*hy
          endif
	    
          x3=(x1+x2)/2.0d0
	      y3=(y1+y2)/2.0d0		! x3,y3是积分边的中点
	    
          perm1=coef(x1,y1) 
          perm2=coef(x2,y2)
          perm3=coef(x3,y3)

          call get_num(n2,node2,i2,j2)
          call get_num(n4,node4,i2,j2)
          do 2000 node=1,4
              bx1=basis_x(node,x1,y1,i1,j1,x0,y0,hx,hy)
              by1=basis_y(node,x1,y1,i1,j1,x0,y0,hx,hy)
              bx2=basis_x(node,x2,y2,i1,j1,x0,y0,hx,hy)
              by2=basis_y(node,x2,y2,i1,j1,x0,y0,hx,hy)
              bx3=basis_x(node,x3,y3,i1,j1,x0,y0,hx,hy)
              by3=basis_y(node,x3,y3,i1,j1,x0,y0,hx,hy)
              value2=value21*(bx1*xn1+by1*yn1)*(1/6.0d0)
     +              +value22*(bx2*xn1+by2*yn1)*(1/6.0d0)
     +              +value23*(bx3*xn1+by3*yn1)*(4/6.0d0)
              value4=value41*(bx1*xn1+by1*yn1)*(1/6.0d0)
     +              +value42*(bx2*xn1+by2*yn1)*(1/6.0d0)
     +              +value43*(bx3*xn1+by3*yn1)*(4/6.0d0)
              call get_num(n1,node,i1,j1)
              
C               print*,'bx1',bx1,'bx2',bx2,'bx3',bx3
C               print*,'by1',by1,'by2',by2,'by3',by3
C               print*,value2,value4
              
              do 1000 m=1,4
              do 1000 l=1,4
              stf(m,l)=stf(m,l)
     +        +p(m,mi1,mj1,n1)*p(l,mi2,mj2,n2)*perm3*value2
     +        +p(m,mi1,mj1,n1)*p(l,mi2,mj2,n4)*perm3*value4
 1000     continue
 2000     continue
 
C           print*,mi1,mj1,i1,j1,nvx,nvy
C           print*,x0,y0,x1,y1,x2,y2,x3,y3
C           do i=1,4
C           do j=1,4
C 		      print*,i,j,stf(i,j)
C           enddo
C           enddo
C           print*,
 
 80         continue

      return
      end
C     -----------------------------------------------
C     FUNCTION penalty term Edge_Integral edge parameter kx=1,ky=0--horizontal;
C        kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse
C        kx1,ky1,kx2,ky2 correponding to edge, equal to nx or ny or 0 
C	have classified
C     -----------------------------------------------
      subroutine p_edge_integ_c(stf,mi1,mj1,kx1,ky1,
     +                mi2,mj2,kx2,ky2,
     +                li1,lj1,nvx,nvy,np,nz,p)
      implicit none
      integer index1,mi1,mj1,kx1,ky1,index2,mi2,mj2,kx2,ky2,
     +         li1,lj1,nvx,nvy,np,nz
      double precision p(4,np,np,nz),stf(4,4)

      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      double precision x1,x2,y1,y2,xyn2,xn1,yn1,perm1,perm2,value1,
     +           value2,basis_x,basis_y,bx,by,perm3
      integer kx,ky,node1,node2,node3,node4,i,j,i1,j1,i2,j2,
     +         k,m,l,n1,n2,n3,n4
     
       stf=0.d0
       xyn2=dsqrt(dble(nvx)**2+dble(nvy)**2)
       xn1=dble(nvx)/xyn2
       yn1=dble(nvy)/xyn2  !normal direction
       kx=abs(nvy)
       ky=abs(nvx)
       if (xn1 .eq. -1 .and. yn1 .eq. 0) then
         node1=1
         node2=1
         node3=4
         node4=4
      else if (xn1 .eq. 1 .and. yn1 .eq. 0) then
         node1=2
         node2=2
	     node3=3
         node4=3
      else if (xn1 .eq. 0 .and. yn1 .eq. -1) then
	     node1=1
         node2=1
	     node3=2
         node4=2
      else if  (xn1 .eq. 0 .and. yn1.eq. 1) then
	     node1=4
         node2=4
	     node3=3
         node4=3
      endif
       do 80 k=1,nx
        i=max(kx*k,1)
        j=max(ky*k,1)
        i1=i+kx1
        j1=j+ky1
        i2=i+kx2
        j2=j+ky2
        value1=hx*xyn2/3.d0
        value2=value1/2.d0
        call get_num(n1,node1,i1,j1)
        call get_num(n2,node2,i2,j2)
        call get_num(n3,node3,i1,j1)
        call get_num(n4,node4,i2,j2)

       do 1000 m=1,4
       do 1000 l=1,4
        stf(m,l)=stf(m,l)
     +    +p(m,mi1,mj1,n1)*p(l,mi2,mj2,n2)*value1
     +    +p(m,mi1,mj1,n1)*p(l,mi2,mj2,n4)*value2
     +    +p(m,mi1,mj1,n3)*p(l,mi2,mj2,n2)*value2
     +    +p(m,mi1,mj1,n3)*p(l,mi2,mj2,n4)*value1
C 		if (mi1.eq.1 .and. mj1.eq.2) then
C 		if (m.eq.1 .and. l.eq. 1) then
C 			print*,p(m,mi1,mj1,n1),p(l,mi2,mj2,n2),value1
C 			print*,p(m,mi1,mj1,n1),p(l,mi2,mj2,n4),value2
C 			print*,p(m,mi1,mj1,n3),p(l,mi2,mj2,n2),value2
C 			print*,p(m,mi1,mj1,n3),p(l,mi2,mj2,n4),value1
C 		endif
C 		endif
1000   continue 
 80    continue
 
C 	  print*,mi1,mj1,nvx,nvy
C 	  print*,node1,node2,node3,node4
C 	  do k=1,4
C 	  do l=1,4
C 		print*,k,l,stf(k,l)
C 	  enddo
C 	  enddo
C 	  print*,
      return
      end
C     -----------------------------------------------
C     FUNCTION Edge_Integral edge parameter  细对粗
C        kx1,ky1,kx2,ky2 correponding to edge, equal to nx-1 or ny-1 or 0 
C	have classified
C     -----------------------------------------------
      subroutine edge_integ_fc(stf,mi1,mj1,ni1,nj1,
     +                mi2,mj2,kx2,ky2,
     +                li1,lj1,nvx,nvy,k,np,nz,p)
      implicit none
      integer index1,mi1,mj1,index2,mi2,mj2,kx2,ky2,
     +         li1,lj1,nvx,nvy,k,np,nz,ni1,nj1
      double precision p(4,np,np,nz),stf(4,4)
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      double precision xyn2,xn1,yn1,perm1,perm2,coef,
     +           basis_x,basis_y,bx1,by1,bx2,by2,bx3,by3,
     +           x0,y0,x1,y1,x2,y2,x3,y3,perm3,
     +           value2,value4,value21,value22,value23,
     +           value41,value42,value43
      integer kx,ky,node2,node4,i,j,i2,j2,ii1,jj1,m,l,
     +         ntran,n2,n4
     
      stf=0.d0
	  
      x0=xleft+dble(mi1-1)*xmax
      y0=yleft+dble(mj1-1)*ymax
	  
      xyn2=dsqrt(dble(nvx)**2+dble(nvy)**2)
      xn1=dble(nvx)/xyn2
      yn1=dble(nvy)/xyn2  !normal direction
      kx=abs(nvy)
      ky=abs(nvx) !kx=1,ky=0--horizontal; kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse
	  
      i=max(kx*k,1)
      j=max(ky*k,1)
      i2=i+kx2
      j2=j+ky2
      ii1=ntran(ni1)
      jj1=ntran(nj1)
	  
      if (xn1 .eq. -1 .and. yn1 .eq. 0) then
	     if (kx2.eq.0 .and. ky2.eq.0) then
	        node2=1
	        node4=4
		    x1=x0+(ii1  )*hx
	        y1=y0+(jj1-1)*hy
	        x2=x0+(ii1  )*hx
	        y2=y0+(jj1  )*hy
	     else
		    node2=2
	        node4=3
	        x1=x0+(ii1-1)*hx
	        y1=y0+(jj1-1)*hy
	        x2=x0+(ii1-1)*hx
	        y2=y0+(jj1  )*hy
	     endif
      else if (xn1 .eq. 1 .and. yn1 .eq. 0) then
	      if (kx2.eq.0 .and. ky2.eq.0) then
             node2=1
	         node4=4
		     x1=x0+(ii1  )*hx
		     y1=y0+(jj1-1)*hy
	         x2=x0+(ii1  )*hx
	         y2=y0+(jj1  )*hy
	      else
		     node2=2
	         node4=3
		     x1=x0+(ii1-1)*hx
		     y1=y0+(jj1-1)*hy
	         x2=x0+(ii1-1)*hx
		     y2=y0+(jj1  )*hy
	      endif
	      
      else if (xn1 .eq. 0 .and. yn1 .eq. -1) then
	      if (kx2.eq.0 .and. ky2.eq.0) then
             node2=1
             node4=2
		     x1=x0+(ii1-1)*hx
		     y1=x0+(jj1  )*hy
		     x2=x0+(ii1  )*hx
		     y2=x0+(jj1  )*hy
	      else
		     node2=4
	         node4=3
		     x1=x0+(ii1-1)*hx
		     y1=y0+(jj1-1)*hy
		     x2=x0+(ii1  )*hx
		     y2=y0+(jj1-1)*hy
	      endif

      else if  (xn1 .eq. 0 .and. yn1.eq. 1) then
	      if (kx2.eq.0 .and. ky2.eq.0) then
             node2=1
             node4=2
		     x1=x0+(ii1-1)*hx
		     y1=y0+(jj1  )*hy
		     x2=x0+(ii1  )*hx
		     y2=y0+(jj1  )*hy
	      else
		     node2=4
	         node4=3
		     x1=x0+(ii1-1)*hx
		     y1=y0+(jj1-1)*hy
		     x2=x0+(ii1  )*hx
		     y2=y0+(jj1-1)*hy
	      endif
      endif
C        print*,kx,ky,node2,node4
	  	  
      x3=(x1+x2)/2.0d0
      y3=(y1+y2)/2.0d0
		  
      perm1=coef(x1,y1)   
      perm2=coef(x2,y2) 
      perm3=coef(x3,y3)

      call get_num(n2,node2,i2,j2)
      call get_num(n4,node4,i2,j2)
      
      value21=hx*(1.0d0  )
      value22=hx*(0.0d0  )
      value23=hx*(1/2.0d0)
      
      value41=hx*(0.0d0  )
      value42=hx*(1.0d0  )
      value43=hx*(1/2.0d0)

      do 1000 m=1,4
	     bx1=basis_x(m,x1,y1,ii1,jj1,x0,y0,hx,hy)
         by1=basis_y(m,x1,y1,ii1,jj1,x0,y0,hx,hy)
         bx2=basis_x(m,x2,y2,ii1,jj1,x0,y0,hx,hy)
         by2=basis_y(m,x2,y2,ii1,jj1,x0,y0,hx,hy)
         bx3=basis_x(m,x3,y3,ii1,jj1,x0,y0,hx,hy)
         by3=basis_y(m,x3,y3,ii1,jj1,x0,y0,hx,hy)
         value2=value21*(bx1*xn1+by1*yn1)*(1/6.0d0)
     +         +value22*(bx2*xn1+by2*yn1)*(1/6.0d0)
     +         +value23*(bx3*xn1+by3*yn1)*(4/6.0d0)
         value4=value41*(bx1*xn1+by1*yn1)*(1/6.0d0)
     +         +value42*(bx2*xn1+by2*yn1)*(1/6.0d0)
     +         +value43*(bx3*xn1+by3*yn1)*(4/6.0d0)
         do 2000 l=1,4
            stf(m,l)=stf(m,l)
     +          +p(l,mi2,mj2,n2)*perm3*value2
     +          +p(l,mi2,mj2,n4)*perm3*value4
     
C  		if (ii1.eq.1 .and. jj1.eq.2) then
C  		if (mi2.eq.1 .and. mj2.eq.2) then
C 		if (mi1.eq.1 .and. mj1.eq.1) then
C 		if (xn1.eq.0 .and. yn1.eq.-1) then
C 			print*,m,l,node2,node4,n2,n4
C 			print*,'bx1',bx1,'by1',by1
C 			print*,'x3',x3,'y3',y3
C  			print*,p(m,mi2,mj2,n2),p(m,mi2,mj2,n4),perm1,value1
C  			print*,p(l,mi2,mj2,n2)*perm1*value1
C      +             +p(l,mi2,mj2,n4)*perm1*value1
C  			print*,
C 		endif
C 		endif
C  		endif
C  		endif
     
 2000     continue
 1000     continue
 
C       print*,mi1,mj1,ni1,nj1,mi2,mj2,nvx,nvy
C       print*,node2,node4,ii1,jj1
C       print*,'x0:',x0,'y0:',y0,'x3:',x3,'y3:',y3
C       print*,'x1:',x1,'y1:',y1,'x2:',x2,'y2:',y2
C       print*,
C       do i=1,4
C       do j=1,4
C 	      print*,i,j,stf(i,j)
C       enddo
C       enddo
 
      return
      end
C     -----------------------------------------------
C     FUNCTION Edge_Integral edge parameter 
C        kx1,ky1,kx2,ky2 correponding to edge, equal to nx-1 or ny-1 or 0 
C	have classified
C     -----------------------------------------------
      subroutine edge_integ_cf(stf,mi1,mj1,kx1,ky1,
     +                mi2,mj2,ni2,nj2,
     +                li1,lj1,nvx,nvy,k,np,nz,p)
      implicit none
      integer index1,mi1,mj1,ni1,nj1,kx1,ky1,index2,mi2,mj2,
     +         ni2,nj2,li1,lj1,nvx,nvy,k,np,nz
      double precision p(4,np,np,nz),stf(4,4)
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft

      double precision xyn2,xn1,yn1,perm1,perm2,coef,perm3,
     +           basis_x,basis_y,bx1,by1,bx2,by2,bx3,by3,
     +           x0,y0,x1,y1,x2,y2,x3,y3,value2,value4,
     +           value21,value41,value22,value42,value23,value43
      integer kx,ky,node2,node4,i,j,ii1,jj1,m,n1,i1,j1,ntran,
     +		  node,l,ii2,jj2
      
        stf=0.d0
	  
	    x0=xleft+dble(mi1-1)*xmax
        y0=yleft+dble(mj1-1)*ymax
	  
        xyn2=dsqrt(dble(nvx)**2+dble(nvy)**2)
        xn1=dble(nvx)/xyn2
        yn1=dble(nvy)/xyn2  !normal direction
        kx=abs(nvy)
        ky=abs(nvx) !kx=1,ky=0--horizontal; kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse
		
        i=max(kx*k,1)
        j=max(ky*k,1)
        i1=i+kx1
        j1=j+ky1
        ii1=ntran(i1)
        jj1=ntran(j1)
		
        if (xn1 .eq. -1 .and. yn1 .eq. 0) then
		    if (kx1.eq.0 .and. ky1.eq.0) then
			    node2=2
		        node4=3
			    x1=x0+(ii1-1)*hx
			    y1=y0+(jj1-1)*hy
			    x2=x0+(ii1-1)*hx
			    y2=y0+(jj1  )*hy
		    else
			    node2=1
			    node4=4
			    x1=x0+(ii1  )*hx
			    y1=y0+(jj1-1)*hy
			    x2=x0+(ii1  )*hx
			    y2=y0+(jj1  )*hy
		    endif
	    
        else if (xn1 .eq. 1 .and. yn1 .eq. 0) then
		    if (kx1.eq.0 .and. ky1.eq.0) then
			    node2=2
			    node4=3
			    x1=x0+(ii1-1)*hx
			    y1=y0+(jj1-1)*hy
			    x2=x0+(ii1-1)*hx
			    y2=y0+(jj1  )*hy
            else
			    node2=1
			    node4=4
			    x1=x0+(ii1  )*hx
			    y1=y0+(jj1-1)*hy
			    x2=x0+(ii1  )*hx
			    y2=y0+(jj1  )*hy
		    endif

        else if (xn1 .eq. 0 .and. yn1 .eq. -1) then
		    if (kx1.eq.0 .and. ky1.eq.0) then
			    node2=4
			    node4=3
			    x1=x0+(ii1-1)*hx
			    y1=y0+(jj1-1)*hy
			    x2=x0+(ii1  )*hx
			    y2=y0+(jj1-1)*hy
		    else
			    node2=1
			    node4=2
			    x1=x0+(ii1-1)*hx
			    y1=y0+(jj1  )*hy
			    x2=x0+(ii1  )*hx
			    y2=y0+(jj1  )*hx
		    endif

        else if  (xn1 .eq. 0 .and. yn1.eq. 1) then
		    if (kx1.eq.0 .and. ky1.eq.0) then
			    node2=4
			    node4=3
			    x1=x0+(ii1-1)*hx
			    y1=y0+(jj1-1)*hy
			    x2=x0+(ii1  )*hx
			    y2=y0+(jj1-1)*hy
		    else
			    node2=1
			    node4=2
			    x1=x0+(ii1-1)*hx
			    y1=y0+(jj1  )*hy
			    x2=x0+(ii1  )*hx
			    y2=y0+(jj1  )*hy
		    endif

        endif
C        print*,kx,ky,node2,node4
        
        x3=(x1+x2)/2.0d0
        y3=(y1+y2)/2.0d0
	  
        perm1=coef(x1,y1)
        perm2=coef(x2,y2)
        perm3=coef(x3,y3)
	  
C         perm1=coef(xleft+(ii1-1)*hx,yleft+(jj1-1)*hy)  !perm(ii1,jj1)
C         perm2=coef(xleft+(ii1+kx-1)*hx,yleft+(jj1+ky-1)*hy) !perm(ii1+kx,jj1+ky)
        value21=hx*(1.0d0  )
        value22=hx*(0.0d0  )
        value23=hx*(1/2.0d0)
        
        value41=hx*(0.0d0  )
        value42=hx*(1.0d0  )
        value43=hx*(1/2.0d0)

C 		print*,'begin'
C 		print*,mi1,mj1,ni1,nj1,mi2,mj2,ni2,nj2,nvx,nvy
        do 2000 node=1,4
	     bx1=basis_x(node,x1,y1,ii1,jj1,x0,y0,hx,hy)
         by1=basis_y(node,x1,y1,ii1,jj1,x0,y0,hx,hy)
         bx2=basis_x(node,x2,y2,ii1,jj1,x0,y0,hx,hy)
         by2=basis_y(node,x2,y2,ii1,jj1,x0,y0,hx,hy)
         bx3=basis_x(node,x3,y3,ii1,jj1,x0,y0,hx,hy)
         by3=basis_y(node,x3,y3,ii1,jj1,x0,y0,hx,hy)
         value2=value21*(bx1*xn1+by1*yn1)*(1/6.0d0)
     +         +value22*(bx2*xn1+by2*yn1)*(1/6.0d0)
     +         +value23*(bx3*xn1+by3*yn1)*(4/6.0d0)
         value4=value41*(bx1*xn1+by1*yn1)*(1/6.0d0)
     +         +value42*(bx2*xn1+by2*yn1)*(1/6.0d0)
     +         +value43*(bx3*xn1+by3*yn1)*(4/6.0d0)
         call get_num(n1,node,i1,j1)
         do 1000 m=1,4
          stf(m,node2)=stf(m,node2)
     +              +p(m,mi1,mj1,n1)*perm3*value2
          stf(m,node4)=stf(m,node4)
     +              +p(m,mi1,mj1,n1)*perm3*value4
		
C 		if (ni2.eq.1 .and. nj2.eq.4) then
C 		if (mi2.eq.2 .and. mj2.eq.1) then
C 			print*,m,node2,node4,x1,y1,x0,y0,xn1,yn1
C 			print*,p(m,mi1,mj1,n1),perm1,value
C 			print*,p(m,mi1,mj1,n1)*perm1*value
C 			print*,
C 		endif
C 		endif
			
 1000    continue
 2000   continue
 
C       print*,mi1,mj1,mi2,mj2,ni2,nj2,nvx,nvy
C       print*,node2,node4,ii1,jj1
C       print*,'x0:',x0,'y0:',y0,'x3:',x3,'y3:',y3
C       print*,'x1:',x1,'y1:',y1,'x2:',x2,'y2:',y2
C       do i=1,4
C       do j=1,4
C 	        print*,i,j,stf(i,j)
C       enddo
C       enddo
C       print*,
      
      return
      end

C     -----------------------------------------------
C     FUNCTION penalty term Edge_Integral 粗对细
C        edge parameter kx=1,ky=0--horizontal;
C        kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse
C        kx1,ky1,kx2,ky2 correponding to edge, equal to nx or ny or 0 
C		 indexi=1,参与计算的是右下两条边
C	     indexi=2,参与计算的是左上两条边
C     -----------------------------------------------
      subroutine p_edge_integ_cf(stf,mi1,mj1,kx1,ky1,
     +        nvx,nvy,k,np,nz,p)
      implicit none
      integer index1,mi1,mj1,kx1,ky1,index2,mi2,mj2,
     +         li1,lj1,nvx,nvy,np,nz
      double precision p(4,np,np,nz),stf(4,4)
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      double precision xyn2,perm1,perm2,value1,perm3,
     +           value2,basis_x,basis_y,bx,by,coef
      integer kx,ky,node1,node2,node3,node4,i,j,i1,j1,
     +        m,n1,n3,k,l
      
      stf=0.d0
c0        x1=xleft+dble(mi1-1)*xmax
c0        y1=yleft+dble(mj1-1)*ymax

        xyn2=dsqrt(dble(nvx)**2+dble(nvy)**2)
        kx=abs(nvy)
        ky=abs(nvx)
        if (nvx .eq. -1 .and. nvy .eq. 0) then
		    if (kx1.eq.0 .and. ky1.eq.0) then
		        node1=1
		        node3=4
				node2=2
				node4=3
	        else
		        node1=2
				node3=3
				node2=1
				node4=4
	        endif
        else if (nvx .eq. 1 .and. nvy .eq. 0) then
	        if (kx1.eq.0 .and. ky1.eq.0) then
			      node1=1
				  node3=4
				  node2=2
				  node4=3
	        else
			      node1=2
				  node3=3
				  node2=1
				  node4=4
	        endif
        else if (nvx .eq. 0 .and. nvy .eq. -1) then
		      if (kx1.eq.0 .and. ky1.eq.0) then
			        node1=1
					node3=2
					node2=4
					node4=3
		      else
			        node1=4
					node3=3
					node2=1
					node4=2
		      endif
        else if  (nvx .eq. 0 .and. nvy.eq. 1) then
	        if (kx1.eq.0 .and. ky1.eq.0) then
			      node1=1
				  node3=2
				  node2=4
				  node4=3
			else
			      node1=4
				  node3=3
				  node2=1
				  node4=2
			endif
        endif
        i=max(kx*k,1)
        j=max(ky*k,1)
        i1=i+kx1
        j1=j+ky1
        value1=hx*xyn2/3.d0
        value2=value1/2.d0
        call get_num(n1,node1,i1,j1)
        call get_num(n3,node3,i1,j1)

C         print*,mi1,mj1,nvx,nvy,i1,j1
C         print*,node2,node4,node1,node3
       do 1000 m=1,4
        stf(m,node2)=stf(m,node2)
     +    +p(m,mi1,mj1,n1)*value1
     +    +p(m,mi1,mj1,n3)*value2
        stf(m,node4)=stf(m,node4)
     +    +p(m,mi1,mj1,n1)*value2
     +    +p(m,mi1,mj1,n3)*value1
1000   continue

C       print*,mi1,mj1,i1,j1
C       do i=1,4
C       do j=1,4
C           print*,i,j,stf(i,j)
C       enddo
C       enddo
C       print*,

      return
      end
C     -----------------------------------------------
C     FUNCTION force right Element_Integral细网格单元上积分(f,\phi)_k
C	have classified
C     -----------------------------------------------
      subroutine r_element_integ_f(fp,mi,mj,ni,nj)
      implicit none
      integer mi,mj,nindex,ni,nj
      double precision fp(4)
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my     
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      double precision volume,force
      integer ii1,jj1,ntran 
      
      volume=hx*hx/4.0d0
      ii1=ntran((mi-1)*nx+ni)
      jj1=ntran((mj-1)*ny+nj)
      fp(1)=force(ii1,jj1)*volume
      fp(2)=force(ii1+1,jj1)*volume
      fp(3)=force(ii1+1,jj1+1)*volume
      fp(4)=force(ii1,jj1+1)*volume
      
      return
      end
C     -----------------------------------------------
C     FUNCTION Edge_Integral edge parameter 
C        kx1,ky1,kx2,ky2 correponding to edge, equal to nx-1 or ny-1 or 0 
C	have classified
C     -----------------------------------------------
      subroutine edge_integ_f(stf,mi,mj,ni,nj,nvx,nvy)
      implicit none
      integer mi,mj,ni,nj,index1,index2,nvx,nvy
      double precision stf(4,4)
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft

      double precision xyn2,xn1,yn1,perm1,perm2,value1,x0,y0,
     +           basis_x,basis_y,bx1,by1,bx2,by2,bx3,by3,
     +           value21,value22,value23,value41,value42,
     +           value43,value2,value4,perm3,
     +           coef,x1,x2,x3,y1,y2,y3
      integer kx,ky,node2,node4,ii1,jj1,m,ntran,k,l
      
      stf=0.d0
	 
      ii1=ni
	  jj1=nj
	    
      x0=xleft+dble(mi-1)*xmax
      y0=yleft+dble(mj-1)*ymax

      xyn2=dsqrt(dble(nvx)**2+dble(nvy)**2)
      xn1=dble(nvx)/xyn2
      yn1=dble(nvy)/xyn2  !normal direction
      kx=abs(nvy)
      ky=abs(nvx) !kx=1,ky=0--horizontal; kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse
      if (xn1 .eq. -1 .and. yn1 .eq. 0) then
	      node2=1
	      node4=4
	      x1=x0+(ii1-1)*hx
	      y1=y0+(jj1-1)*hy
	      x2=x0+(ii1-1)*hx
	      y2=y0+(jj1  )*hy
      else if (xn1 .eq. 1 .and. yn1 .eq. 0) then
          node2=2
	      node4=3
	      x1=x0+(ii1  )*hx
	      y1=y0+(jj1-1)*hy
	      x2=x0+(ii1  )*hx
	      y2=y0+(jj1  )*hy
      else if (xn1 .eq. 0 .and. yn1 .eq. -1) then
          node2=1
          node4=2
          x1=x0+(ii1-1)*hx
	      y1=y0+(jj1-1)*hy
	      x2=x0+(ii1  )*hx
	      y2=y0+(jj1-1)*hy
      else if  (xn1 .eq. 0 .and. yn1.eq. 1) then
          node2=4
          node4=3
          x1=x0+(ii1-1)*hx
	      y1=y0+(jj1  )*hy
	      x2=x0+(ii1  )*hx
	      y2=y0+(jj1  )*hy
      endif
	  
      perm1=coef(x1,y1)
      perm2=coef(x2,y2)
      x3=(x1+x2)/2.0d0
      y3=(y1+y2)/2.0d0
      perm3=coef(x3,y3)
	  
C       value=hx*xyn2/2.d0
      value21=hx*xyn2*(1.0d0  )
      value23=hx*xyn2*(1/2.0d0)
      value22=hx*xyn2*(0.0d0  )
      value41=hx*xyn2*(0.0d0  )
      value43=hx*xyn2*(1/2.0d0)
      value42=hx*xyn2*(1.0d0  )

C 	  print*,mi,mj,ni,nj,nvx,nvy,ii1,jj1
C 	  print*,x0,y0,x1,y1,x2,y2,x3,y3
      
      do 1000 m=1,4
      
	     bx1=basis_x(m,x1,y1,ii1,jj1,x0,y0,hx,hy)
         by1=basis_y(m,x1,y1,ii1,jj1,x0,y0,hx,hy)
         bx2=basis_x(m,x2,y2,ii1,jj1,x0,y0,hx,hy)
         by2=basis_y(m,x2,y2,ii1,jj1,x0,y0,hx,hy)
         bx3=basis_x(m,x3,y3,ii1,jj1,x0,y0,hx,hy)
         by3=basis_y(m,x3,y3,ii1,jj1,x0,y0,hx,hy)
         value2=value21*(bx1*xn1+by1*yn1)*(1/6.0d0)
     +         +value22*(bx2*xn1+by2*yn1)*(1/6.0d0)
     +         +value23*(bx3*xn1+by3*yn1)*(4/6.0d0)
         value4=value41*(bx1*xn1+by1*yn1)*(1/6.0d0)
     +         +value42*(bx2*xn1+by2*yn1)*(1/6.0d0)
     +         +value43*(bx3*xn1+by3*yn1)*(4/6.0d0)
         stf(m,node2)=perm3*value2
         stf(m,node4)=perm3*value4
       
 1000     continue
       
C        print*,mi,mj,ni,nj,nvx,nvy
C        print*,node2,node4
C        print*,x0,y0,x1,y1,x2,y2,x3,y3
C        do k=1,4
C        do l=1,4
C            print*,k,l,stf(k,l)
C        enddo
C        enddo
C        print*,
       
      return
      end

C     -----------------------------------------------
C     FUNCTION penalty term Edge_Integral 细网格单元上加罚项计算 J^0(u,v)
C         edge parameter kx=1,ky=0--horizontal;
C        kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse
C        kx1,ky1,kx2,ky2 correponding to edge, equal to nx or ny or 0 
C	have clasified
C     -----------------------------------------------
      subroutine p_edge_integ_f(stf,nvx,nvy)
      implicit none
      integer index1,index2,nvx,nvy
      double precision stf(4,4)
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my,k,l
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      double precision xyn2,xn1,yn1,value1,value2
      integer kx,ky,node1,node2,node3,node4
      stf=0.d0
      xyn2=dsqrt(dble(nvx)**2+dble(nvy)**2)
      xn1=dble(nvx)/xyn2
      yn1=dble(nvy)/xyn2  !normal direction
      kx=abs(nvy)
      ky=abs(nvx)
      if (xn1 .eq. -1 .and. yn1 .eq. 0) then
         node1=1
         node2=1
         node3=4
         node4=4
      else if (xn1 .eq. 1 .and. yn1 .eq. 0) then
         node1=2
         node2=2
	     node3=3
         node4=3
      else if (xn1 .eq. 0 .and. yn1 .eq. -1) then
	     node1=1
         node2=1
	     node3=2
         node4=2
      else if  (xn1 .eq. 0 .and. yn1.eq. 1) then
	     node1=4
         node2=4
	     node3=3
         node4=3
      endif
      value1=hx*xyn2/3.d0
      value2=value1/2.d0

      stf(node1,node2)=value1
      stf(node1,node4)=value2
      stf(node3,node2)=value2
      stf(node3,node4)=value1
		 
C       print*,nvx,nvy
C       do k=1,4
C       do l=1,4
C          print*,k,l,stf(k,l)
C       enddo
C       enddo
C       print*,

      return
      end
      
      
C   -------------------------------------------------
C   J1 penalty
C
C   -------------------------------------------------
    
    
C     -----------------------------------------------
C     FUNCTION J1 penalty term Edge_Integral edge parameter kx=1,ky=0--horizontal;
C        kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse
C        kx1,ky1,kx2,ky2 correponding to edge, equal to nx or ny or 0 
C     -----------------------------------------------     
      subroutine pj1_edge_integ_c(stf,mi1,mj1,kx1,ky1,
     +                mi2,mj2,kx2,ky2,
     +                li1,lj1,nvx,nvy,np,nz,p)
      implicit none
      integer index1,mi1,mj1,kx1,ky1,index2,mi2,mj2,kx2,ky2,
     +         li1,lj1,nvx,nvy,np,nz
      double precision p(4,np,np,nz),stf(4,4)
      integer nxs,nxend 
       
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      
      double precision xyn2,xn1,yn1,perm1,perm2,value1,value2,
     +           value,basis_x,basis_y,bx,by,coef,perm3,
     +           x,y,x0,y0,coor_gauss
      integer kx,ky,node1,i,j,i1,j1,i2,j2,ii1,jj1,
     +         ntran,m,l,k,n1,n2,node
     
        stf=0.d0

        x0 = xleft + dble(mi1-1)*xmax
        y0 = yleft + dble(mj1-1)*ymax

        xyn2=dsqrt(dble(nvx)**2+dble(nvy)**2)
        xn1=dble(nvx)/xyn2  !单位法向量
        yn1=dble(nvy)/xyn2  !normal direction
        kx=abs(nvy) !kx=1,ky=0--水平; kx=0,ky=1,垂直; kx=1,ky=1,斜边
        ky=abs(nvx) !kx=1,ky=0--horizontal; kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse

       do 80 k=1,nx
        i=max(kx*k,1)
        j=max(ky*k,1)
        i1=i+kx1
        j1=j+ky1
        i2=i+kx2
        j2=j+ky2
        ii1=ntran((li1-1)*nx+i)
        jj1=ntran((lj1-1)*ny+j)
        
        x=x0+(i1-0.5d0)*hx
        y=y0+(j1-0.5d0)*hy
        
        perm1=coef(x,y)**2

        value=hx*xyn2*perm1 !细网格上相邻两点(ii1,jj1) (ii1+kx,jj1+ky)之间的积分，积分区间长度 hx*xyn2

        
        do 2000 node=1,4
         x=coor_gauss(node,1,i1,j1,x0,y0)
         y=coor_gauss(node,2,i1,j1,x0,y0)
         bx=basis_x(node,x,y,i1,j1,x0,y0,hx,hy)
         by=basis_y(node,x,y,i1,j1,x0,y0,hx,hy)
         value1=value*(bx*xn1+by*yn1)
         call get_num(n1,node,i1,j1)
         do 2000 node1=1,4
           bx=basis_x(node1,x,y,i1,j1,x0,y0,hx,hy)
           by=basis_y(node1,x,y,i1,j1,x0,y0,hx,hy)
           value2=value1*(bx*xn1+by*yn1)
           call get_num(n2,node1,i2,j2)
           do 1000 m=1,4
           do 1000 l=1,4        
            stf(m,l)=stf(m,l)
     +        +p(m,mi1,mj1,n1)*p(l,mi2,mj2,n2)*value2!/4.0d0
c0   +    +p(index1,m,mi1,mj1,n1)*p2(l,node2)*perm1*value
c0   +    +p(index1,m,mi1,mj1,n1)*p2(l,node4)*perm2*value
 1000     continue
 2000     continue
 
 80         continue
 
        stf=0.0d0
        
      return
      end      
      
      
C     -----------------------------------------------
C     FUNCTION J1 penalty term Edge_Integral 粗对细
C        kx1,ky1,kx2,ky2 correponding to edge, equal to nx-1 or ny-1 or 0 
C     -----------------------------------------------
      subroutine pj1_edge_integ_cf(stf,index1,mi1,mj1,kx1,ky1,
     +                index2,mi2,mj2,li1,lj1,
     +                nvx,nvy,k,np,nz,p)
      implicit none
      integer index1,mi1,mj1,kx1,ky1,index2,mi2,mj2,
     +         li1,lj1,nvx,nvy,k,np,nz
      double precision p(4,np,np,nz),stf(4,4)
      integer nxs,nxend
              
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      double precision xyn2,xn1,yn1,perm1,perm2,perm3,value,value1,
     +           value2,basis_x,basis_y,bx,by,coef,x01,y01,
     +           x1,y1,x2,y2,x02,y02,coor_gauss,x,y
      integer kx,ky,node,node1,i,j,ii1,jj1,m,n1,i1,j1,ntran
      
        stf=0.d0
        
        x01=xleft+dble(mi1-1)*xmax
        y01=yleft+dble(mj1-1)*ymax
        x02=xleft+dble(mi2-1)*xmax
        y02=yleft+dble(mj2-1)*ymax
        
        xyn2=dsqrt(dble(nvx)**2+dble(nvy)**2)
        xn1=dble(nvx)/xyn2
        yn1=dble(nvy)/xyn2  !normal direction
        kx=abs(nvy)
        ky=abs(nvx) !kx=1,ky=0--horizontal; kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse

        i=max(kx*k,1)
        j=max(ky*k,1)
        i1=i+kx1
        j1=j+ky1
        ii1=ntran((li1-1)*nx+i)
        jj1=ntran((lj1-1)*ny+j)
        
        x1=xleft+(ii1-0.5d0)*hx
        y1=yleft+(jj1-0.5d0)*hy
        x2=xleft+(ii1+kx-0.5d0)*hx
        y2=yleft+(jj1+ky-0.5d0)*hy
        
        perm1=coef(x1,y1)**2
        perm2=coef(x2,y2)**2
        
        value=hx*xyn2*perm1*perm2 !细网格上相邻两点(ii1,jj1) (ii1+kx,jj1+ky)之间的积分，积分区间长度 hx*xyn2

        do 2000 node=1,4
         x=coor_gauss(node,1,i,j,x01,y01)
         y=coor_gauss(node,2,i,j,x01,y01)
         bx=basis_x(node,x,y,i,j,x01,y01,hx,hy)
         by=basis_y(node,x,y,i,j,x01,y01,hx,hy)
         value1=value*(bx*xn1+by*yn1)
         call get_num(n1,node,i1,j1)
         do 2000 node1=1,4
           bx=basis_x(node1,x2,y2,ii1+kx,jj1+ky,x02,y02,hx,hy)
           by=basis_y(node1,x2,y2,ii1+kx,jj1+ky,x02,y02,hx,hy)
           value2=value1*(bx*xn1+by*yn1)
           do 1000 m=1,4    
            stf(m,node1)=stf(m,node1)
     +        +p(m,mi1,mj1,n1)*value2
 1000     continue
 2000     continue
 
         stf=0.0d0
        
      return
      end
      
      
C     -----------------------------------------------
C     FUNCTION penalty term Edge_Integral 细网格单元上加罚项计算 J^1(u,v)
C         edge parameter kx=1,ky=0--horizontal;
C        kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse
C        kx1,ky1,kx2,ky2 correponding to edge, equal to nx or ny or 0 
C     -----------------------------------------------      
      subroutine pj1_edge_integ_f(stf,mi,mj,ni,nj,
     +             nvx,nvy)
      implicit none
      integer mi,mj,ni,nj,nvx,nvy
      double precision stf(4,4)
      integer nxs,nxend
                  
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      double precision xyn2,xn1,yn1,perm1,perm2,perm3,value,value1,
     +           value2,basis_x,basis_y,bx,by,coef,x,y,x0,y0,coor_gauss
      integer kx,ky,node,node1,ii1,jj1,m,ntran,i1,j1
      
        stf=0.d0
        
        ii1=ni
        jj1=nj
        
        x0=xleft+dble(mi-1)*xmax
        y0=yleft+dble(mj-1)*ymax
        
        xyn2=dsqrt(dble(nvx)**2+dble(nvy)**2)
        xn1=dble(nvx)/xyn2
        yn1=dble(nvy)/xyn2  !normal direction
        kx=abs(nvy)
        ky=abs(nvx) !kx=1,ky=0--horizontal; kx=0,ky=1,vertical; kx=1,ky=1,hypotenuse

        i1=ntran((mi-1)*nx+ni)
        j1=ntran((mj-1)*ny+nj)
        
        x=xleft+(i1-0.5d0)*hx
        y=yleft+(j1-0.5d0)*hy
        
        perm1=coef(x,y)**2
        
        value=hx*xyn2*perm1

        do 2000 node=1,4
          x=coor_gauss(node,1,ni,nj,x0,y0)
          y=coor_gauss(node,2,ni,nj,x0,y0)
          bx=basis_x(node,x,y,ii1,jj1,x0,y0,hx,hy)
          by=basis_y(node,x,y,ii1,jj1,x0,y0,hx,hy)
          value1=value*(bx*xn1+by*yn1)!*coef(x,y)     
          do 2000 node1=1,4
            bx=basis_x(node1,x,y,ii1,jj1,x0,y0,hx,hy)
            by=basis_y(node1,x,y,ii1,jj1,x0,y0,hx,hy)
            value2=value1*(bx*xn1+by*yn1)!*coef(x,y)
            stf(node,node1)=stf(node,node1)+value2
 2000     continue

        !stf=0.0d0
        
      return
      end
