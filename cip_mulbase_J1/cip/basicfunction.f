C   
C   2015年12月31日,修订的刚度矩阵装配函数
C   可以很好的解决常系数问题，变系数问题有待改善
C
C   此程序包包含粗细网格上基函数的计算，细网格上单刚计算，
C   单元编号的计算，自由度编号，局部整体编号之间映射关系，边界自由度赋值
C   节点坐标计算，边界条件 r，扩散系数定义 coef， 源项 f
C****************************************************
C     计算单元编号
C     返回 k---当前单元（mi,mj)（+(ni,nj)）的单元编号
C****************************************************
      subroutine getel_num(k,mxmax,i_fnum,mxy,mi,mj,
     &                    nindex,ni,nj)
      implicit none
      integer i_fnum(mxy),mxmax(my),k,mxy,mindex,mi,mj,nindex,ni,nj   
        
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/ nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      integer k_fine,k_coarse,k0,j

      k_fine=0  ! 单元 (mi,mj) 之前所有细网格单元的总数
      k_coarse=0 ! 单元 (mi,mj) 之前所有粗网格单元的总数
      k0=0         !宏单元的编号
      do j=1,mj-1
        k0=k0+mxmax(j)
      enddo
      k0=k0+mi-1  ! (mi,mj) 宏单元的编号-1
C     print*,k0,i_fnum(k0),nx,mx,mi,mj
c     k0=(mj-1)*mx+mi-1
      if(k0.gt.0) then
        k_fine=nx*ny*i_fnum(k0)   ! before (mi,mj), total number of elements of fine mesh
        k_coarse=k0-i_fnum(k0)   !before (mi,mj), total number of elements of coarse mesh
      endif
      if(nindex.eq.0) then
        k=k_coarse+k_fine+1 !当前网格为粗网格
      else
        k=k_coarse+k_fine+(nj-1)*nx+ni !当前网格为细网格
C          ((nj-1)*nx+ni-1)*2+3-nindex-1 计算的是 (ni,nj) 之前的单元数   
      endif 
      return
      end
C********************************************************************      
C     assemble the dof, produce the local-global dof map
C     按自由度编号顺序，产生编号映射关系   
C********************************************************************
      subroutine assembledof(lgmap,nel,mxmax,i_fnum,mxy,iflag,
     &            kdof,kbdof,kidof,ub,kbmax)
      implicit none
      integer i_fnum(mxy),mxmax(my),lgmap(nel,4),
     &   iflag(-1:mx+1,-1:my+1),nel,mxy,kdof,kbdof,kidof,kbmax
      double precision ub(kbmax)
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/ nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      integer mi,mj,i,j,ni,nj,ktemp,kel,itrue,imove
      double precision x0,y0,r,basis_x,basis_y
      
      kdof=0
      kidof=kbmax   !内部节点自由度编号
      kbdof=0   !边界自由度编号
      do 20 mj=1,my
       do 10 mi=1,mxmax(mj)
        x0=xleft+(mi-1)*xmax
        y0=yleft+(mj-1)*ymax
        if(iflag(mi,mj).eq.0) then !粗网格单元
          kdof=kdof+1
          call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,0,0,0)
          if(iflag(mi-1,mj).lt.0.or.iflag(mi,mj-1).lt.0) then
            kbdof=kbdof+1
            ub(kbdof)=r(x0,y0)
            ktemp=kbdof
          else
            kidof=kidof+1
            ktemp=kidof
          endif
          lgmap(kel,1)=ktemp
C1                    
          itrue=0
	    
C	    左下单元
          if(iflag(mi-1,mj-1).eq.0) then
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj-1,0,0,0)
            lgmap(kel,3)=ktemp
          elseif(iflag(mi-1,mj-1).eq.1)then		!加密单元
            itrue=1
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj-1,2,nx,ny)
            lgmap(kel,3)=ktemp+1
          endif    
	    
C	     正下方单元
          if(iflag(mi,mj-1).eq.0) then
            call getel_num(kel,mxmax,i_fnum,mxy,mi,mj-1,0,0,0)
            lgmap(kel,4)=ktemp
          elseif(iflag(mi,mj-1).eq.1) then
            itrue=1
            call getel_num(kel,mxmax,i_fnum,mxy,mi,mj-1,2,1,ny)
            lgmap(kel,4)=ktemp+1
          endif
	    
C	    左单元
          if(iflag(mi-1,mj).eq.0) then
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj,0,0,0)
            lgmap(kel,2)=ktemp
          elseif(iflag(mi-1,mj).eq.1) then
            itrue=1
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj,1,nx,1)
            lgmap(kel,2)=ktemp+1
          endif

C	    在上述几个单元中出现了加密单元
          if(itrue.eq.1) then
	       kdof=kdof+1
             if(iflag(mi-1,mj).lt.0.or.iflag(mi,mj-1).lt.0) then
                 kbdof=kbdof+1
                 ub(kbdof)=r(x0,y0)
             else
                 kidof=kidof+1
             endif
          endif 
	    
        else  !加密网格单元
          nj=1
          ni=1
          kdof=kdof+1
          call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,2,ni,nj)
          if(iflag(mi-1,mj).lt.0.or.iflag(mi,mj-1).lt.0) then
             kbdof=kbdof+1
             ub(kbdof)=r(x0,y0)             
             ktemp=kbdof
          else
             kidof=kidof+1
             ktemp=kidof
          endif
          lgmap(kel,1)=ktemp
C2                   
          itrue=0
	    
C		左下单元
          if(iflag(mi-1,mj-1).eq.0) then
            itrue=1
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj-1,0,0,0)
            lgmap(kel,3)=ktemp+1
          elseif(iflag(mi-1,mj-1).eq.1) then
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj-1,2,nx,ny)
            lgmap(kel,3)=ktemp
          endif
	    
C		下方单元
          if(iflag(mi,mj-1).eq.0) then
            itrue=1
            call getel_num(kel,mxmax,i_fnum,mxy,mi,mj-1,0,0,0)
            lgmap(kel,4)=ktemp+1
          elseif(iflag(mi,mj-1).eq.1) then
            call getel_num(kel,mxmax,i_fnum,mxy,mi,mj-1,2,1,ny)
            lgmap(kel,4)=ktemp
          endif
	    
C		左方单元
          if(iflag(mi-1,mj).eq.0) then
            itrue=1
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj,0,0,0)
            lgmap(kel,2)=ktemp+1
          elseif(iflag(mi-1,mj).eq.1) then
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj,1,nx,1)
            lgmap(kel,2)=ktemp

          endif
	    
	    
          if(itrue.eq.1) then
            kdof=kdof+1
            if(iflag(mi-1,mj).lt.0.or.iflag(mi,mj-1).lt.0) then
              kbdof=kbdof+1
              ub(kbdof)=r(x0,y0)              
            else
              kidof=kidof+1
            endif
          endif          
C2
          nj=1
          do ni=2,nx
            kdof=kdof+1
            call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,2,ni,nj)
            if(iflag(mi,mj-1).lt.0) then
              kbdof=kbdof+1
              ub(kbdof)=r(x0+(ni-1)*hx,y0)              
              ktemp=kbdof
            else
              kidof=kidof+1
              ktemp=kidof
            endif
            lgmap(kel-1,2)=ktemp
            lgmap(kel,1)=ktemp
            if(iflag(mi,mj-1).eq.1) then
              call getel_num(kel,mxmax,i_fnum,mxy,mi,mj-1,2,ni,ny)
              lgmap(kel-1,3)=ktemp
              lgmap(kel,4)=ktemp
            endif
          enddo
C2          

C		右边单元未加密
          if(iflag(mi+1,mj).lt.1) then
            ni=nx
            do nj=1,ny-1
              kdof=kdof+1
              call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,2,ni,nj)
              if(iflag(mi+1,mj).lt.0) then
                kbdof=kbdof+1
                ub(kbdof)=r(x0+hx*nx,y0+nj*hy)              
                ktemp=kbdof
              else
                kidof=kidof+1
                ktemp=kidof
              endif
              lgmap(kel,3)=ktemp
              lgmap(kel+nx,2)=ktemp
            enddo
          endif
C2          
          if(iflag(mi,mj+1).lt.1) then
            nj=ny
            do ni=1,nx-1
              kdof=kdof+1
              call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,1,ni,nj)
              if(iflag(mi,mj+1).lt.0) then
                kbdof=kbdof+1
                ub(kbdof)=r(x0+hx*ni,y0+hy*ny)              
                ktemp=kbdof
              else
                kidof=kidof+1
                ktemp=kidof
              endif
              lgmap(kel,3)=ktemp
              lgmap(kel+1,4)=ktemp
            enddo
          endif
C2         
          ni=1
          do nj=2,ny
            kdof=kdof+1
            call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,2,ni,nj)
            if(iflag(mi-1,mj).lt.0) then
              kbdof=kbdof+1
              ub(kbdof)=r(x0,y0+hy*(nj-1))              
              ktemp=kbdof
            else
              kidof=kidof+1
              ktemp=kidof
            endif            
            lgmap(kel,1)=ktemp
            lgmap(kel-nx,4)=ktemp
            if(iflag(mi-1,mj).eq.1) then
              call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj,1,nx,nj)
              lgmap(kel,2)=ktemp
              lgmap(kel-nx,3)=ktemp
            endif
          enddo
C2          
          do nj=2,ny
          do ni=2,nx
            kdof=kdof+1
            kidof=kidof+1
            call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,2,ni,nj)
            lgmap(kel,1)=kidof
            lgmap(kel-nx,4)=kidof
            lgmap(kel-1,2)=kidof
            lgmap(kel-nx-1,3)=kidof         
          enddo
          enddo         
        endif
10     continue
C      处理右端边界粗网格点
       mi=mxmax(mj)
       if(iflag(mi,mj).eq.0) then !粗网格单元
         kdof=kdof+1
         kbdof=kbdof+1
         ub(kbdof)=r(x0+xmax,y0)        
         call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,0,0,0)
         lgmap(kel,2)=kbdof
         itrue=0
         if(iflag(mi,mj-1).eq.0) then
           call getel_num(kel,mxmax,i_fnum,mxy,mi,mj-1,0,0,0)
           lgmap(kel,3)=kbdof
         elseif(iflag(mi,mj-1).eq.1) then
           itrue=1
           call getel_num(kel,mxmax,i_fnum,mxy,mi,mj-1,2,nx,ny)
           lgmap(kel,3)=kbdof+1
         endif
         if(itrue.eq.1)then
           kdof=kdof+1
           kbdof=kbdof+1
           ub(kbdof)=r(x0+xmax,y0)
         endif        
       else  !加密网格单元
         kdof=kdof+1
         kbdof=kbdof+1
         ub(kbdof)=r(x0+xmax,y0)        
         call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,1,nx,1)
         lgmap(kel,2)=kbdof
         itrue=0
         if(iflag(mi,mj-1).eq.0) then
           itrue=1
           call getel_num(kel,mxmax,i_fnum,mxy,mi,mj-1,0,0,0)
           lgmap(kel,3)=kbdof+1
         elseif(iflag(mi,mj-1).eq.1) then
           call getel_num(kel,mxmax,i_fnum,mxy,mi,mj-1,2,nx,ny)
           lgmap(kel,3)=kbdof
         endif
         if(itrue.eq.1)then
           kdof=kdof+1
           kbdof=kbdof+1
           ub(kbdof)=r(x0+xmax,y0)
         endif
       endif
20    continue
C  处理上边界粗网格点
      mj=my
      do mi=1,mxmax(mj)
        x0=xleft+(mi-1)*xmax
        y0=yleft+(mj-1)*ymax
        if(iflag(mi,mj).eq.0) then !粗网格单元
          kdof=kdof+1
          kbdof=kbdof+1
          ub(kbdof)=r(x0,y0+ymax)
          call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,0,0,0)
          lgmap(kel,4)=kbdof
          itrue=0
          if(iflag(mi-1,mj).eq.0) then
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj,0,0,0)
            lgmap(kel,3)=kbdof
          elseif(iflag(mi-1,mj).eq.1) then
            itrue=1
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj,2,nx,ny)
            lgmap(kel,3)=kbdof+1
          endif
          if(itrue.eq.1)then
            kdof=kdof+1
            kbdof=kbdof+1
            ub(kbdof)=r(x0,y0+ymax)
          endif
        else  !加密网格单元
          kdof=kdof+1
          kbdof=kbdof+1
          ub(kbdof)=r(x0,y0+ymax)
          call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,2,1,ny)
          lgmap(kel,4)=kbdof
          itrue=0
          if(iflag(mi-1,mj).eq.0) then
            itrue=1
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj,0,0,0)
            lgmap(kel,3)=kbdof+1
          elseif(iflag(mi-1,mj).eq.1) then
            call getel_num(kel,mxmax,i_fnum,mxy,mi-1,mj,2,nx,ny)
            lgmap(kel,3)=kbdof
          endif
          if(itrue.eq.1)then
            kdof=kdof+1
            kbdof=kbdof+1
            ub(kbdof)=r(x0,y0+ymax)
          endif
        endif
      enddo
C  处理右上角网格点
      mj=my
      mi=mxmax(mj)
      x0=xleft+(mi-1)*xmax
      y0=yleft+(mj-1)*ymax
      if(iflag(mi,mj).eq.0) then !粗网格单元
        kdof=kdof+1
        kbdof=kbdof+1
        ub(kbdof)=r(x0+xmax,y0+ymax)
        call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,0,0,0)
        lgmap(kel,3)=kbdof
      else  !加密网格单元
        kdof=kdof+1
        kbdof=kbdof+1
        ub(kbdof)=r(x0+xmax,y0+ymax)        
        call getel_num(kel,mxmax,i_fnum,mxy,mi,mj,2,nx,ny)
        lgmap(kel,3)=kbdof
      endif      
        
C    对内部节点重新编号，移动 kbmax-kbdof

      imove=kbmax-kbdof
      kidof=kidof-kbmax
      do i=1,nel
      do j=1,4
        if(lgmap(i,j).gt.kbdof)lgmap(i,j)=lgmap(i,j)-imove
      enddo
      enddo
C    检验编号
      if(kbdof+kidof.ne.kdof)print*,'Wrong DOF'  
C      print*,kbdof,kidof,kdof,imove    
	
	open(1,file='lg',status='unknown')
	write(1,*) ((lgmap(i,j),j=1,4),i=1,nel)
	close(1)
    
      return
      end
C     -------------------------------------------------
C     SUBROUTINE GET_NUM 获得粗网格单元内细网格节点的编号
C     --------------------------------------------------
      subroutine get_num(num,node,i,j)
      implicit none
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
C     ------------------------------------------------
C     FUNCTION NTRAN
C     ------------------------------------------------
      function ntran(i)
      implicit none
      integer i,ntran 
      
      ntran=i 
      return
      end
C    -----------------------------------------------
C     SUBROUTINE GET_STIFF 细网格上单刚
C     -----------------------------------------------
      subroutine get_stiff(stf,mi,mj,i,j)
      implicit none
      double precision stf(4,4)
      integer index,mi,mj,i,j
       
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my
      common /param/nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft

      double precision x,y,permc,coef,bk1,bk2,bl1,bl2,area,
     +                 basis_x,basis_y,x0,y0,coor_gauss,d,
     +                 xtmp,ytmp
      integer k,l,node


      area=hx*hy
      x0=xleft+dble(mi-1)*xmax
      y0=yleft+dble(mj-1)*ymax
	
      xtmp=x0+dble(i-0.5)*hx
      ytmp=y0+dble(j-0.5)*hy
	
C       stf(1,1)=(2.0d0/3.0d0)*coef(xtmp,ytmp)
C       stf(1,2)=-(1.0d0/6.0d0)*coef(xtmp,ytmp)
C       stf(1,3)=-(1.0d0/3.0d0)*coef(xtmp,ytmp)
C       stf(1,4)=-(1.0d0/6.0d0)*coef(xtmp,ytmp)
C       stf(2,1)=-(1.0d0/6.0d0)*coef(xtmp,ytmp)
C       stf(2,2)=(2.0d0/3.0d0)*coef(xtmp,ytmp)
C       stf(2,3)=-(1.0d0/6.0d0)*coef(xtmp,ytmp)
C       stf(2,4)=-(1.0d0/3.0d0)*coef(xtmp,ytmp)
C       stf(3,1)=-(1.0d0/3.0d0)*coef(xtmp,ytmp)
C       stf(3,2)=-(1.0d0/6.0d0)*coef(xtmp,ytmp)
C       stf(3,3)=(2.0d0/3.0d0)*coef(xtmp,ytmp)
C       stf(3,4)=-(1.0d0/6.0d0)*coef(xtmp,ytmp)
C       stf(4,1)=-(1.0d0/6.0d0)*coef(xtmp,ytmp)
C       stf(4,2)=-(1.0d0/3.0d0)*coef(xtmp,ytmp)
C       stf(4,3)=-(1.0d0/6.0d0)*coef(xtmp,ytmp)
C       stf(4,4)=(2.0d0/3.0d0)*coef(xtmp,ytmp)

      do k=1,4
         do l=1,4
            stf(k,l)=0.0d0
            do node=1,4
               x=coor_gauss(node,1,i,j,x0,y0)
               y=coor_gauss(node,2,i,j,x0,y0)
               bk1=basis_x(k,x,y,i,j,x0,y0,hx,hy)
               bk2=basis_y(k,x,y,i,j,x0,y0,hx,hy)
               bl1=basis_x(l,x,y,i,j,x0,y0,hx,hy)
               bl2=basis_y(l,x,y,i,j,x0,y0,hx,hy)
               stf(k,l)=stf(k,l)+(area/4.0d0)*
     +              (d(1,1,x,y)*bk1*bl1+d(1,2,x,y)*bk1*bl2
     +              +d(2,1,x,y)*bk2*bl1+d(2,2,x,y)*bk2*bl2)
            end do
         end do
      end do

      return
      end

	

C     --------------------------------------------
C     FUNCTION COOR
C     --------------------------------------------
      real*8 function coor_gauss(node,kcoor,i,j,xl,yl)
      implicit none
      double precision coor_gauss,x0,y0,xl,yl
      integer node,kcoor,i,j
       
      double precision hx,hy
      integer nx,ny  
      common /param/ nx,ny,hx,hy

      real*8 tmp,xleft,yleft,xmax,ymax

      x0=xl+(i-1)*hx
      y0=yl+(j-1)*hy
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
C     FUNCTION BASIS_X
C     ----------------------------------------------
      real*8 function basis_x(node,x,y,i,j,xl,yl,hx,hy)
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      real*8 hx,hy,x,y,x0,y0,area,xl,yl,xmax,ymax
      integer nx,ny,index,node,kcoor,i,j

      x0=xl+(i-1)*hx
      y0=yl+(j-1)*hy
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
      real*8 function basis_y(node,x,y,i,j,xl,yl,hx,hy)
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      real*8 hx,hy,x,y,x0,y0,area,xl,yl,xmax,ymax
      integer nx,ny,index,node,kcoor,i,j

      x0=xl+(i-1)*hx
      y0=yl+(j-1)*hy
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

C ********************************************************
C     SOURCE FUNCTION
C********************************************************
      function force(i,j)
      implicit none
      double precision force
      integer i,j
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my      
      common /param/ nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      double precision x,y
            
      x=xleft+dble(i-1)*hx
      y=yleft+dble(j-1)*hy
      force=1.0d0
c      force=-10.0d0
c      force=1.0d0
      end

C********************************************************
C     EXACT SOLUTION for boundary
C ********************************************************
      function r(x,y)
      implicit none
      double precision r,x,y
      
      double precision xmax,ymax,hx,hy,xleft,yleft
      integer nx,ny,mx,my 
      common /param/ nx,ny,hx,hy
      common /domain/ mx,my,xmax,ymax,xleft,yleft
      
      double precision pi
            
      pi=dacos(-1.d0)
C      x=xleft+dble(i-1)*hx
C      y=yleft+dble(j-1)*hy
c      if(x.lt.0.d0)then
c        r=(x**2+y**2)**(1.d0/3.d0)*dsin(2.d0/3.d0*(datan(y/x)+pi))
c      elseif(x.gt.0.d0.and.y.lt.0.d0) then
c        r=(x**2+y**2)**(1.d0/3.d0)*dsin(2.d0/3.d0*(datan(y/x)+pi*2.d0))
c      elseif(x.ge.0.d0.and.y.gt.0) then
c        r=(x**2+y**2)**(1.d0/3.d0)*dsin(2.d0/3.d0*datan(y/x))
c      else
c        r=0.d0
c      endif
c        r=x**2+y**2
c        r=(1.d0-x)*x*(1.d0-y)*y/2.d0!+2.d0
         r=0.d0
      return
      end
C********************************************************
C     读取多尺度基函数
C********************************************************      
      subroutine openbase(p,mx,my,nz)
      implicit none
      double precision p(4,mx,my,nz)
      integer mx,my,nz,mi,mj,nxyf,i
      
      open(1,file='dat/base/nbase1',status='old')
      open(2,file='dat/base/nbase2',status='old')
      open(3,file='dat/base/nbase3',status='old')
      open(4,file='dat/base/nbase4',status='old')

      do 150 mj=1,my 
      do 150 mi=1,mx
C.....READ THE BASIS INFORMATION
      read(1,*) nxyf, (p(1,mi,mj,i),i=1,nxyf)
      read(2,*) nxyf, (p(2,mi,mj,i),i=1,nxyf)
      read(3,*) nxyf, (p(3,mi,mj,i),i=1,nxyf) 
      read(4,*) nxyf, (p(4,mi,mj,i),i=1,nxyf)

150   continue

      close(1)
      close(2)
      close(3)
      close(4)
      end             
C********************************************************
C  刚度矩阵元素值合并，仍然存放在原来矩阵中
C*********************************************************
      subroutine arrangeam(head,m2,value)    
      use typedef
      implicit none
      type(row_el), pointer :: head,p,q
      integer m2,itrue,irr
      double precision value
      
      p=>head
      q=>head    
      itrue=0

      do while (associated(p))
        if(p%col.eq.m2) then
          p%value=p%value+value
          itrue=1
          exit
        else
C          print*,'1',associated(p%next)
          q=>p
          p=>p%next
        endif
      end do

      if(itrue.eq.0) then
C        print*,'2'
        allocate(p,stat=irr)
        if(irr.ne.0) print*,'Wrong allocate in arrange'
        p%col=m2
        p%value=value
        p%next=>null()
        q%next=>p
      endif
      
      return
      end
      
C***************************************************
C    计算刚度矩阵 按照单元顺序计算，保存在稀疏矩阵 a 中
C***************************************************
      subroutine assemblestiff(a,drhs,kdof,lgmap,nel,mxmax,i_fnum,mxy,
     &            np,nz,p,iflag,kbdof,delta2,symm)
      use typedef
      implicit none
      type(pointerArray),pointer :: a(:)
      integer lgmap(nel,4),mxmax(my),i_fnum(mxy)
      integer np,nz,kdof,kbdof,nel,mxy,iflag(-1:mx+1,-1:my+1)
      double precision delta2,delta3,drhs(kdof),symm,p(4,np,np,nz)
      
      double precision xmax,ymax,hx,hy,xleft,yleft,eps
      integer nx,ny,mx,my
      common /param/ nx,ny,hx,hy   !细网格尺寸，参数
      common /domain/ mx,my,xmax,ymax,xleft,yleft !粗网格尺寸，参数
      common/epsilon/eps  !周期系数参数
      
      integer mi,mj,m,m1,m2,l,k1,ni,nj,nindex,k2,k
      double precision stf1(4,4),stf2(4,4),stf3(4,4),stf4(4,4),value,
     +                   r_element_integ_c,fp1(4),fp2(4)
      
      interface
        subroutine  arrangeam(am,m2,value)
          use typedef
          implicit none
          type(row_el),pointer :: am
          integer m2
          double precision value
        end subroutine
      end interface
      
      delta3=0.1d0*hx
      
      do 20 mj=1,my
      do 10 mi=1,mxmax(mj)
C        print*,'mi,mj=',mi,mj
        if(iflag(mi,mj).eq.0) then !粗网格单元
          call getel_num(k1,mxmax,i_fnum,mxy,mi,mj,0,0,0)
          call element_integ_c(stf1,mi,mj,mx,nz,p)
          do m=1,4
            m1=lgmap(k1,m)
C             if (m1.eq.16) then
C                 print*,1,k1
C                 do k=1,4
C                     print*,lgmap(k1,k)
C                     print*,m,k,stf1(m,k)
C                 enddo
C                 print*,
C             endif
            if(m1.gt.kbdof) then
              do l=1,4
                m2=lgmap(k1,l)
                value=stf1(m,l)
                call arrangeam(a(m1)%el,m2,value)
              enddo
            endif
          enddo
          
C         form the right hand side
          do m=1,4
            drhs(lgmap(k1,m))=drhs(lgmap(k1,m))
     +          +r_element_integ_c(mi,mj,m,mx,nz,p)
          enddo
	    
C2        下面考虑粗细网格单元交界边界积分 
          if(iflag(mi-1,mj).eq.1) then  ! left
C           粗网格单元内
C 		    stf1=0.0d0
C 		    stf2=0.0d0
C 		    stf3=0.0d0
            call edge_integ_c(stf1,mi,mj,0,0,mi,mj,0,0,
     +               mi,mj,-1,0,mx,nz,p)
            call p_edge_integ_c(stf2,mi,mj,0,0,mi,mj,0,0,
     +               mi,mj,-1,0,mx,nz,p)
            call pj1_edge_integ_c(stf3,mi,mj,0,0,mi,mj,0,0,
     +               mi,mj,-1,0,mx,nz,p)
            do m=1,4
              m1=lgmap(k1,m)
C               if (m1.eq.16) then
C                 print*,2,k1
C                 print*,m,'symm:',symm,'delta2:',delta2
C                 do l=1,4
C                     print*,lgmap(k1,l)
C                     print*,l,stf1(l,m),stf1(m,l),stf2(l,m)
C                 enddo
C                 print*,
C               endif
              if(m1.gt.kbdof) then
                do l=1,4
                  m2=lgmap(k1,l)
                  value=-(stf1(m,l)+symm*stf1(l,m))/2.d0
     +                 +delta2*stf2(l,m)+delta3*stf3(l,m)        !加上stf3(l,m)
                  call arrangeam(a(m1)%el,m2,value)
                enddo
              endif
            enddo

C           粗网格单元与左边细网格单元
            ni=nx
            do 101 nj=1,ny
C 		      stf1=0.0d0
C 		      stf2=0.0d0
C 		      stf3=0.0d0
              call edge_integ_fc(stf1,mi-1,mj,ni,nj,mi,mj,0,0, 
     +               mi,mj,-1,0,nj,mx,nz,p)
              !print*,'begin',1
              call edge_integ_cf(stf2,mi,mj,0,0,mi-1,mj,ni,nj,
     +               mi,mj,-1,0,nj,mx,nz,p)
              !print*,'end',1
              call p_edge_integ_cf(stf3,mi,mj,0,0,-1,0,nj,
     &               mx,nz,p)
              call pj1_edge_integ_cf(stf4,2,mi,mj,0,0,1,mi-1,mj,
     +               mi,mj,-1,0,nj,mx,nz,p)
              call getel_num(k2,mxmax,i_fnum,mxy,mi-1,mj,
     &               1,ni,nj)
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,3,'k1:',k1,'k2:',k2
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k2,l)
C                         print*,l,stf1(l,m),stf2(m,l),stf3(m,l)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4
                    m2=lgmap(k2,l)
C                     value=-(stf1(l,m)+symm*stf2(m,l))/2.d0
C      +                  -delta2*stf3(m,l)
                    value=-(-stf2(m,l)+symm*stf1(l,m))/2.0d0
     -                    -delta2*stf3(m,l)-delta3*stf4(m,l)                    
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
 101        continue            
          endif 
          
          if(iflag(mi,mj-1).eq.1) then  ! bottom
C           粗网格单元内    
C 		    stf1=0.0d0
C 		    stf2=0.0d0
C 		    stf3=0.0d0
            call edge_integ_c(stf1,mi,mj,0,0,mi,mj,0,0,
     +                mi,mj,0,-1,mx,nz,p)
            call p_edge_integ_c(stf2,mi,mj,0,0,mi,mj,0,0,
     +                mi,mj,0,-1,mx,nz,p)
            call pj1_edge_integ_c(stf3,mi,mj,0,0,mi,mj,0,0,
     +                mi,mj,0,-1,mx,nz,p)
            call getel_num(k1,mxmax,i_fnum,mxy,mi,mj,
     &                0,0,0)              
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,4,'k1:',k1
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k1,l)
C                         print*,l,stf1(l,m),stf1(m,l),stf2(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4             
                    m2=lgmap(k1,l)
                    value=-(stf1(m,l)+symm*stf1(l,m))/2.d0 
     +                 +delta2*stf2(l,m)+delta3*stf3(l,m)             !加上stf3(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
C           粗网格单元与下边细网格单元              
            nj=ny
            do 102 ni=1,nx
C 		    stf1=0.0d0
C 		    stf2=0.0d0
C 		    stf3=0.0d0
              call edge_integ_fc(stf1,mi,mj-1,ni,nj,mi,mj,0,0,
     +                mi,mj,0,-1,ni,mx,nz,p)
              !print*,'begin',2
              call edge_integ_cf(stf2,mi,mj,0,0,mi,mj-1,ni,nj,
     +                mi,mj,0,-1,ni,mx,nz,p)
              !print*,'end',2
              call p_edge_integ_cf(stf3,mi,mj,0,0,0,-1,ni,
     &               mx,nz,p)
              call pj1_edge_integ_cf(stf4,1,mi,mj,0,0,2,mi,mj-1,
     +                mi,mj,0,-1,ni,mx,nz,p)
              call getel_num(k2,mxmax,i_fnum,mxy,mi,mj-1,
     &           1,ni,nj)
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,5,'k1:',k1,'k2:',k2
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k2,l)
C                         print*,l,stf1(l,m),stf2(m,l),stf3(m,l)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4             
                    m2=lgmap(k2,l)
C                     value=-(stf1(l,m)+symm*stf2(m,l))/2.d0
C      +                 -delta2*stf3(m,l)
                    value=-(-stf2(m,l)+symm*stf1(l,m))/2.0d0
     -                    -delta2*stf3(m,l)-delta3*stf4(m,l)   
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
 102        continue
          endif               

          if(iflag(mi+1,mj).eq.1) then  ! right
C           粗网格单元内 
C 		    stf1=0.0d0
C 		    stf2=0.0d0
C 		    stf3=0.0d0
            call edge_integ_c(stf1,mi,mj,nx-1,0,mi,mj,nx-1,0,
     +         mi+1,mj,1,0,mx,nz,p)
            call p_edge_integ_c(stf2,mi,mj,nx-1,0,mi,mj,nx-1,0,
     +         mi+1,mj,1,0,mx,nz,p)
            call pj1_edge_integ_c(stf3,mi,mj,nx-1,0,mi,mj,nx-1,0,
     +         mi+1,mj,1,0,mx,nz,p)
            call getel_num(k1,mxmax,i_fnum,mxy,mi,mj,
     &         0,0,0)              
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,6,'k1:',k1
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k1,l)
C                         print*,l,stf1(l,m),stf1(m,l),stf2(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4             
                    m2=lgmap(k1,l)
                    value=-(stf1(m,l)+symm*stf1(l,m))/2.d0 
     +                   +delta2*stf2(l,m)+delta3*stf3(l,m)       !加上stf3(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
C           粗网格单元与右边细网格单元       
            ni=1
            do 103 nj=1,ny
C 		     stf1=0.0d0
C 		     stf2=0.0d0
C 		     stf3=0.0d0
              call edge_integ_fc(stf1,mi+1,mj,ni,nj,mi,mj,nx-1,0,
     +         mi+1,mj,1,0,nj,mx,nz,p)
              !print*,'begin',3
              call edge_integ_cf(stf2,mi,mj,nx-1,0,mi+1,mj,ni,nj,
     +         mi+1,mj,1,0,nj,mx,nz,p)
              !print*,'end',3
              call p_edge_integ_cf(stf3,mi,mj,nx-1,0,1,0,nj,
     &               mx,nz,p)
              call pj1_edge_integ_cf(stf4,1,mi,mj,nx-1,0,2,mi+1,mj,
     +         mi+1,mj,1,0,nj,mx,nz,p)
              call getel_num(k2,mxmax,i_fnum,mxy,mi+1,mj,
     &           1,ni,nj)
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,7,'k1:',k1,'k2:',k2
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k2,l)
C                         print*,l,stf1(l,m),stf2(m,l),stf3(m,l)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4             
                    m2=lgmap(k2,l)
C                      value=-(stf1(l,m)+symm*stf2(m,l))/2.d0 
C      -                 -delta2*stf3(m,l)
                    value=-(-stf2(m,l)+symm*stf1(l,m))/2.0d0
     -                    -delta2*stf3(m,l)-delta3*stf4(m,l)   
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
103         continue 
          endif          
               
          if(iflag(mi,mj+1).eq.1) then  ! upper
C           粗网格单元内   
C 		    stf1=0.0d0
C 		    stf2=0.0d0
C 		    stf3=0.0d0
            call edge_integ_c(stf1,mi,mj,0,ny-1,mi,mj,0,ny-1,
     +                mi,mj+1,0,1,mx,nz,p)
            call p_edge_integ_c(stf2,mi,mj,0,ny-1,mi,mj,0,ny-1,
     +                mi,mj+1,0,1,mx,nz,p)
            call pj1_edge_integ_c(stf3,mi,mj,0,ny-1,mi,mj,0,ny-1,
     +                mi,mj+1,0,1,mx,nz,p)
            call getel_num(k1,mxmax,i_fnum,mxy,mi,mj,
     &                0,0,0)              
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,8,'k1:',k1
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k1,l)
C                         print*,l,stf1(l,m),stf1(m,l),stf2(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4             
                    m2=lgmap(k1,l)
                    value=-(stf1(m,l)+symm*stf1(l,m))/2.d0
     +                 +delta2*stf2(l,m)+delta3*stf3(l,m)   !加上stf3(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
		  
C           粗网格单元与上边细网格单元       
            nj=1
            do 104 ni=1,nx
C 		    stf1=0.0d0
C 		    stf2=0.0d0
C 		    stf3=0.0d0
              call edge_integ_fc(stf1,mi,mj+1,ni,nj,mi,mj,0,ny-1,
     +                mi,mj+1,0,1,ni,mx,nz,p)
              !print*,'begin',4
              call edge_integ_cf(stf2,mi,mj,0,ny-1,mi,mj+1,ni,nj,
     +                mi,mj+1,0,1,ni,mx,nz,p)
              !print*,'end',4
              call p_edge_integ_cf(stf3,mi,mj,0,ny-1,0,1,ni,
     &               mx,nz,p)
              call pj1_edge_integ_cf(stf4,2,mi,mj,0,ny-1,1,mi,mj+1,
     +                mi,mj+1,0,1,ni,mx,nz,p)
              call getel_num(k2,mxmax,i_fnum,mxy,mi,mj+1,
     &               1,ni,nj)
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,9,'k1:',k1,'k2:',k2
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k2,l)
C                         print*,l,stf1(l,m),stf2(m,l),stf3(m,l)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4             
                    m2=lgmap(k2,l)
C                     value=-(stf1(l,m)+symm*stf2(m,l))/2.d0
C      -                      -delta2*stf3(m,l)
                    value=-(-stf2(m,l)+symm*stf1(l,m))/2.0d0
     -                    -delta2*stf3(m,l)-delta3*stf4(m,l)   
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
 104        continue
          endif                       
C1
        else  !加密网格单元
          do nj=1,ny
          do ni=1,nx
C 		    stf1=0.0d0
C 		    stf2=0.0d0
C 		    stf3=0.0d0
            call getel_num(k1,mxmax,i_fnum,mxy,mi,mj,2,ni,nj)
            call get_stiff(stf1,mi,mj,ni,nj)
            do m=1,4
              m1=lgmap(k1,m)
C               if (m1.eq.16) then
C                     print*,10,'k1:',k1
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k1,l)
C                         print*,l,stf1(m,l)
C                     enddo
C                     print*,
C                 endif
              if(m1.gt.kbdof) then
                do l=1,4
                  m2=lgmap(k1,l)
                  value=stf1(m,l)
                  call arrangeam(a(m1)%el,m2,value)
                enddo
              endif
            enddo            
C           form the right hand side
            call r_element_integ_f(fp1,mi,mj,ni,nj)
            do m=1,4
              drhs(lgmap(k1,m))=drhs(lgmap(k1,m))+fp1(m)
            enddo
          enddo
          enddo
C2        下面考虑粗细网格单元交界边界积分      
          if(iflag(mi-1,mj).eq.0) then  ! left
            call getel_num(k2,mxmax,i_fnum,mxy,mi-1,mj,0,0,0)     
            ni=1
            do 201 nj=1,ny
C             细网格单元内  
C 		      stf1=0.0d0
C 		      stf2=0.0d0
C 		      stf3=0.0d0
              call edge_integ_f(stf1,mi,mj,ni,nj,-1,0)  
              call p_edge_integ_f(stf2,-1,0)
              call pj1_edge_integ_f(stf3,mi,mj,ni,nj,-1,0)
              call getel_num(k1,mxmax,i_fnum,mxy,mi,mj,1,ni,nj)
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,11,'k1:',k1
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k1,l)
C                         print*,l,stf1(l,m),stf1(m,l),stf2(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4
                    m2=lgmap(k1,l)
                    value=-(stf1(m,l)+symm*stf1(l,m))/2.d0
     +                 +delta2*stf2(l,m)+delta3*stf3(l,m)              !加上stf3(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo

C             细网格单元与左边粗网格单元
C 		      stf1=0.0d0
C 		      stf2=0.0d0
C 		      stf3=0.0d0
              !print*,'begin',5
              call edge_integ_cf(stf1,mi-1,mj,nx-1,0,mi,mj,ni,nj,
     +               mi,mj,-1,0,nj,mx,nz,p)
              !print*,'end',5
              call edge_integ_fc(stf2,mi,mj,ni,nj,mi-1,mj,nx-1,0,
     +               mi,mj,-1,0,nj,mx,nz,p)
              call p_edge_integ_cf(stf3,mi-1,mj,nx-1,0,-1,0,nj,
     &               mx,nz,p)
              call pj1_edge_integ_cf(stf4,1,mi-1,mj,nx-1,0,2,mi,mj, 
     +               mi,mj,-1,0,nj,mx,nz,p)
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,12,'k1:',k1,'k2:',k2
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k2,l)
C                         print*,l,stf1(l,m),stf2(m,l),stf3(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4
                    m2=lgmap(k2,l)
C                     value=-(stf1(l,m)+symm*stf2(m,l))/2.d0
C      -               -delta2*stf3(l,m)
                    value=-(-stf2(m,l)+symm*stf1(l,m))/2.d0
     -                 -delta2*stf3(l,m)-delta3*stf4(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
 201        continue
          endif 
          
          if(iflag(mi,mj-1).eq.0) then  ! bottom
            call getel_num(k2,mxmax,i_fnum,mxy,mi,mj-1,0,0,0)
            nj=1
            do 202 ni=1,nx
C             细网格单元内        
C 		      stf1=0.0d0
C 		      stf2=0.0d0
C 		      stf3=0.0d0
              call edge_integ_f(stf1,mi,mj,ni,nj,0,-1)
              call p_edge_integ_f(stf2,0,-1)
              call pj1_edge_integ_f(stf3,mi,mj,ni,nj,0,-1)
              call getel_num(k1,mxmax,i_fnum,mxy,mi,mj,
     &               1,ni,nj)              
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,13,'k1:',k1
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k1,l)
C                         print*,l,stf1(l,m),stf1(m,l),stf2(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4
                    m2=lgmap(k1,l)
                    value=-(stf1(m,l)+symm*stf1(l,m))/2.d0
     +                 +delta2*stf2(l,m)+delta3*stf3(l,m)               !加上stf3(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo

C             细网格单元与下边粗网格单元 
C 		      stf1=0.0d0
C 		      stf2=0.0d0
C 		      stf3=0.0d0
              !print*,'begin',6
              call edge_integ_cf(stf1,mi,mj-1,0,ny-1,mi,mj,ni,nj,
     +                mi,mj,0,-1,ni,mx,nz,p)
              !print*,'end',6
              call edge_integ_fc(stf2,mi,mj,ni,nj,mi,mj-1,0,ny-1,
     +                mi,mj,0,-1,ni,mx,nz,p)
              call p_edge_integ_cf(stf3,mi,mj-1,0,ny-1,0,-1,ni,
     &                mx,nz,p)
              call pj1_edge_integ_cf(stf4,2,mi,mj-1,0,ny-1,1,mi,mj,
     +                mi,mj,0,-1,ni,mx,nz,p)
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,14,'k1:',k1,'k2:',k2
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k2,l)
C                         print*,l,stf1(l,m),stf2(m,l),stf3(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4
                    m2=lgmap(k2,l)
C                     value=-(stf1(l,m)+symm*stf2(m,l))/2.d0
C      -               -delta2*stf3(l,m)
                    value=-(-stf2(m,l)+symm*stf1(l,m))/2.d0
     -                 -delta2*stf3(l,m)-delta3*stf4(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
 202        continue
          endif               

          if(iflag(mi+1,mj).eq.0) then  ! right
            call getel_num(k2,mxmax,i_fnum,mxy,mi+1,mj,0,0,0)
            ni=nx
            do 203 nj=1,ny
C             细网格单元内 
C 		      stf1=0.0d0
C 		      stf2=0.0d0
C 		      stf3=0.0d0
              call edge_integ_f(stf1,mi,mj,ni,nj,1,0)
              call p_edge_integ_f(stf2,1,0)
              call pj1_edge_integ_f(stf3,mi,mj,ni+1,nj,1,0)
              call getel_num(k1,mxmax,i_fnum,mxy,mi,mj,
     &              1,ni,nj)              
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,15,'k1:',k1
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k1,l)
C                         print*,l,stf1(l,m),stf1(m,l),stf2(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4
                    m2=lgmap(k1,l)
                    value=-(stf1(m,l)+symm*stf1(l,m))/2.d0
     +                 +delta2*stf2(l,m)+delta3*stf3(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo

C             细网格单元与右边粗网格单元
C 	          stf1=0.0d0
C 		      stf2=0.0d0
C 		      stf3=0.0d0
              call edge_integ_cf(stf1,mi+1,mj,0,0,mi,mj,ni,nj,
     +           mi+1,mj,1,0,nj,mx,nz,p)
              call edge_integ_fc(stf2,mi,mj,ni,nj,mi+1,mj,0,0,
     +           mi+1,mj,1,0,nj,mx,nz,p)
              call p_edge_integ_cf(stf3,mi+1,mj,0,0,1,0,nj,
     &           mx,nz,p)
              call pj1_edge_integ_cf(stf4,2,mi+1,mj,0,0,1,mi,mj,
     +           mi+1,mj,1,0,nj,mx,nz,p)
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,16,'k1:',k1,'k2:',k2
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k2,l)
C                         print*,l,stf1(l,m),stf2(m,l),stf3(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4
                    m2=lgmap(k2,l)
C                     value=-(stf1(l,m)+symm*stf2(m,l))/2.d0
C      -                 -delta2*stf3(l,m)
                    value=-(-stf2(m,l)+symm*stf1(l,m))/2.d0
     -                 -delta2*stf3(l,m)-delta3*stf4(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
 203        continue 
          endif          
               
          if(iflag(mi,mj+1).eq.0) then  ! upper
            call getel_num(k2,mxmax,i_fnum,mxy,mi,mj+1,0,0,0)
            nj=ny
            do 204 ni=1,nx
C             细网格单元内
C 		      stf1=0.0d0
C 		      stf2=0.0d0
              call edge_integ_f(stf1,mi,mj,ni,nj,0,1)
              call p_edge_integ_f(stf2,0,1)
              call pj1_edge_integ_f(stf3,mi,mj,ni,nj+1,0,1)
              call getel_num(k1,mxmax,i_fnum,mxy,mi,mj,
     &                   2,ni,nj)              
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,17,'k1:',k1
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k1,l)
C                         print*,l,stf1(l,m),stf1(m,l),stf2(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4
                    m2=lgmap(k1,l)
                    value=-(stf1(m,l)+symm*stf1(l,m))/2.d0
     +              +delta2*stf2(l,m)+delta3*stf3(l,m)               !加上stf3(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo

C             细网格单元与上边粗网格单元
C 		      stf1=0.0d0
C 		      stf2=0.0d0
C 		      stf3=0.0d0
              !print*,'begin',8
              call edge_integ_cf(stf1,mi,mj+1,0,0,mi,mj,ni,nj,
     +                mi,mj+1,0,1,ni,mx,nz,p)
              call edge_integ_fc(stf2,mi,mj,ni,nj,mi,mj+1,0,0,
     +                mi,mj+1,0,1,ni,mx,nz,p)
              call p_edge_integ_cf(stf3,mi,mj+1,0,0,0,1,ni,
     &                mx,nz,p)
              call pj1_edge_integ_cf(stf4,1,mi,mj+1,0,0,2,mi,mj,
     +                mi,mj+1,0,1,ni,mx,nz,p)
              !print*,'end',8
              do m=1,4
                m1=lgmap(k1,m)
C                 if (m1.eq.16) then
C                     print*,18,'k1:',k1,'k2:',k2
C                     print*,m,'symm:',symm,'delta2:',delta2
C                     do l=1,4
C                         print*,lgmap(k2,l)
C                         print*,l,stf1(l,m),stf2(m,l),stf3(l,m)
C                     enddo
C                     print*,
C                 endif
                if(m1.gt.kbdof) then
                  do l=1,4
                    m2=lgmap(k2,l)
C                     value=-(stf1(l,m)+symm*stf2(m,l))/2.d0
C      -                 -delta2*stf3(l,m)
                    value=-(-stf2(m,l)+symm*stf1(l,m))/2.d0
     -                 -delta2*stf3(l,m)-delta3*stf4(l,m)
                    call arrangeam(a(m1)%el,m2,value)
                  enddo
                endif
              enddo
 204        continue
          endif               
        endif      
10    continue
20    continue
      end    