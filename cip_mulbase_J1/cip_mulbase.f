C   2015年11月09日，四边形剖分的MsFEM
C   粗细网格非协调剖分---带悬点的CIPFEM，线性元   
C   注意粗细网格惩罚系数的匹配  程序名中最后的 c--粗 f-- 细
C***************************************************
      module typedef
        implicit none
        type :: row_el   ! 结构体
          integer :: col
          double precision :: value
          type(row_el),pointer :: next
        end type row_el
        type pointerArray   !指针数组，指针类型结构体
          type(row_el),pointer::el
        end type pointerArray
      end module typedef
C---------------------------------------------------
      program main
      use typedef
      implicit none

C      dimension mxmax(max_ny)  ! L型区域参数 mxmax 每一列最大行标
c      dimension da(3,max_nm*12),drhs(max_nm)  !刚度矩阵，右端项
c      dimension x(max_nm),Ax(max_nm*12)
C      dimension iflag(max_nx,max_ny),i_fnum(max_nx*max_ny) !粗细网格标识
C      integer NAi(max_nm*12), NAp(max_nm+1), index1(max_nm) !调用 UMFPACK 参数
C      动态分配数组
      double precision, allocatable :: p(:,:,:,:),ub(:),u(:)
      double precision, allocatable :: x(:),Ax(:),da(:),db(:),drhs(:)
      integer, allocatable :: mxmax(:),ida(:),jda(:),idb(:),jdb(:)
      integer, allocatable :: iflag(:,:),i_fnum(:),lgmap(:,:)
      integer, allocatable :: NAi(:),NAp(:),index1(:)
C---------------------------------------------------
      type(pointerArray),pointer :: a(:)   !指针数组，指针类型结构体
      type(row_el),pointer :: head,ptmp
C---------------------------------------------------      
      double precision sigma1,sigma2,delta0,delta1,delta2,
     +                  delta3,symm,rho
      integer lmax,ifl,mxy,ineps    
      integer nel,kbmax,limax,lbmax,k1,k2,i,j,k,nz   
      integer mi,mj,kmaxi,kmaxb,kbdof,kidof,kdof
	integer*8 kidof8,limax8
C---------------------------------------------------      
      double precision xmax,ymax,hx,hy,xleft,yleft,xright,yright,eps
      integer nx,ny,mx,my,np     
      common /param/ nx,ny,hx,hy   !细网格尺寸，参数
      common /domain/ mx,my,xmax,ymax,xleft,yleft !粗网格尺寸，参数
      common/epsilon/eps  !周期系数参数
C---------------------------------------------------          
      interface      
        subroutine assemblestiff(a,drhs,kdof,lgmap,nel,mxmax,
     &                i_fnum,mxy,np,nz,p,iflag,kbdof,delta2,symm)
        use typedef
          implicit none
          type(pointerArray),pointer :: a(:)
          integer lgmap(nel,4),mxmax(my),i_fnum(mxy)
          integer np,nz,kdof,kbdof,nel,mxy,iflag(-1:mx+1,-1:my+1)
          double precision delta2,drhs(kdof),symm,p(4,nx,ny,nz)
          double precision xmax,ymax,hx,hy,xleft,yleft
          integer nx,ny,mx,my
          common /param/ nx,ny,hx,hy
          common /domain/ mx,my,xmax,ymax,xleft,yleft
        end subroutine
      end interface
C---------------------------------------------------
C.....PARAMETERS..............................
      open(1,file='inputlr.dat',status='old')
      read(1,*) nx,mx,lmax,ifl,sigma1,sigma2,symm
      read(1,*)
      read(1,*) xleft,yleft,xright,yright
      close(1)
      open (10,file='eps.dat',status='Old')
      read(10,*) ineps
      close(10)      
      print*,nx,mx,xleft,yleft,xright,yright,lmax,ifl,sigma1,sigma2
C     16,64,2,1,100.d0,100.d0,-1.d0,64
C     nx--粗网格内加密个数  mx--粗网格数
C     lmax 奇点附近局部加密单元的层数, ifl=1,代表有局部加密，ifl=0，无局部加密
C     sigma1,sigma2 for penalty,ineps=1/eps
C     symm=-1.0 represents symmetric
C     0.d0,0.d0,1.d0,1.d0	! xleft,yleft,xright,yright----  domain
      ny=nx
      my=mx
      allocate(mxmax(my))
      allocate(iflag(-1:mx+1,-1:my+1))
      mxy=mx*my
      allocate(i_fnum(mxy))
C---------------------------------------------------
      xmax=(xright-xleft)/mx  !粗网格尺寸
      ymax=(yright-yleft)/my
      hx=xmax/nx  !细网格尺寸
      hy=ymax/ny
      eps=1.d0/ineps
      rho=eps
      print*,'first',eps
      delta0=sigma1/xmax  !粗网格加罚系数
      delta1=sigma2/hx    !细网格加罚系数
      delta2=sigma1/rho       !粗细网格交叉时加罚系数
      delta3=sigma2*rho
C---------------------------------------------------
C     Define the shape 可以定义L型区域
      do j=1,my/2
        mxmax(j)=mx   !/2
      enddo
      do j=my/2+1,my
        mxmax(j)=mx
      enddo
C---------------------------------------------------
      do j=-1,my+1
      do i=-1,mx+1
        iflag(i,j)=-1  !刻画区域外的单元 -1
      enddo
      enddo
      do j=1,my
      do i=1,mxmax(j)
        iflag(i,j)=0  !0 代表未加密， 1 代表加密
      enddo
      enddo
	
      do i=1,lmax
        do j=1,my
           iflag(i,j)=ifl
           iflag(j,i)=ifl
	       iflag(mx+1-i,j)=ifl
           iflag(j,my+1-i)=ifl           
        enddo
      enddo
C       iflag(2,2)=ifl

C--------------------------------------------      
      k=0
      do j=1,my
        k1=0
        do i=1,j-1
          k1=k1+mxmax(i) !记录第一层到第j-1层的宏单元数（粗网格）
        enddo
        do i=1,mxmax(j)
          k=k+iflag(i,j)
          k2=k1+i       !宏单元编号，先行后列，即先水平后垂直方向
          i_fnum(k2)=k  !保存的是到第KK1个宏单元之前已加密的宏单元个数
        enddo
      enddo
C---------------------------------------------------      
C     计算总的单元数 nel （粗细网格总数）      
      if(iflag(mxmax(my),my).eq.0) then !判断右上角单元是否加密
        call getel_num(nel,mxmax,i_fnum,mxy,mxmax(my),my,0,0,0)
      else
        call getel_num(nel,mxmax,i_fnum,mxy,mxmax(my),
     &                   my,1,nx,ny)
      endif
	print*,nel
      kbmax=mx*nx*4   !边界自由度编号最大值    4条边
      allocate(lgmap(nel,4))
      allocate(ub(kbmax))
      ub=0.d0
      lgmap=0
C     给每个单元赋自由度编号，边界内部节点分开来编号，同时给边界自由度赋值（Dirichlet）
      call assembledof(lgmap,nel,mxmax,i_fnum,mxy,iflag,
     &        kdof,kbdof,kidof,ub,kbmax)
    
      allocate(drhs(kdof))  !右端项
      allocate(u(kdof))     !解
      drhs=0.d0
      u=0.d0
C---------------------------------------------------      
      nz=(nx+1)*(ny+1)
      allocate(p(4,mx,my,nz))     
      call openbase(p,mx,my,nz)
C---------------------------------------------------      
      allocate(a(1:kdof))
      do k=1,kdof
        allocate(a(k)%el)
        a(k)%el%col=k
        a(k)%el%value=0.d0
        a(k)%el%next=>null()
      end do     
C---------------------------------------------------
c      irefine=1
c      do while (irefine.ne.0)
C---------------------------------------------------
C    生成总刚，右端；将刚度矩阵每一个非0元素保存在a中，右端 drhs
       call assemblestiff(a,drhs,kdof,lgmap,nel,mxmax,i_fnum,mxy,
     +      mx,nz,p,iflag,kbdof,delta2,symm)
c       call assemblestiff(a,drhs,kdof,lgmap,nel,mxmax,i_fnum,mxy,
c     &                mx,nz,p,iflag,kbdof,delta2,delta3,symm)
C---------------------------------------------------
C     将刚度矩阵每一个非0元素保存在da中，ida行号；jda列号；da值
C     将刚度矩阵转换为内部结点存储da 和边界结点存储 db
      kmaxi=kidof*12
      kmaxb=kbdof*12
      allocate(da(kmaxi))   !刚度矩阵  
      allocate(ida(kmaxi))     
      allocate(jda(kmaxi))
      allocate(db(kmaxb))
      allocate(idb(kmaxb))
      allocate(jdb(kmaxb))
      limax=0
      lbmax=0
      do 100 k1=kbdof+1,kdof
        head=>a(k1)%el
        do while (associated(head))
          k2=head%col
          if(k2.gt.kbdof) then
            limax=limax+1
            ida(limax)=k1-kbdof
            jda(limax)=k2-kbdof
            da(limax)=head%value
          else
            lbmax=lbmax+1
            idb(lbmax)=k1-kbdof
            jdb(lbmax)=k2
            db(lbmax)=head%value
          end if
          ptmp=>head
          head=>head%next
          deallocate(ptmp)
        end do
C       deallocate(a(k1)%el)
 100  continue  
      deallocate(a)
C---------------------------------------------------      
      do k=1,kidof
        drhs(k)=drhs(k+kbdof)
      enddo
      do k=1,lbmax
        k1=idb(k)
        k2=jdb(k)
        drhs(k1)=drhs(k1)-db(k)*ub(k2)
      enddo
      deallocate(idb)
      deallocate(jdb)
      deallocate(db)
C---------------------------------------------------      
C     输出刚度矩阵，右端
      open(2,file='../dat/da_cip.dat')
      write(2,*)limax,kidof,kbdof
      do k=1,limax
        write(2,*) ida(k),jda(k),da(k)
      enddo
      do k=1,kbdof
        write(2,*) k,ub(k)
      enddo
      do k=1,kidof
        write(2,*) k,drhs(k)
      enddo
      close(2)
C---------------------------------------------------      
      allocate(x(kidof))
      allocate(Ax(kidof*12))
      allocate(NAi(kidof*12))
      allocate(NAp(kidof+1))
      allocate(index1(kidof))
C---------------------------------------------------      
      print*,nel,kidof,kbdof,kdof,limax,lbmax
	kidof8=kidof
	limax8=limax
      print*, 'go to Umfpack'
      call umfdglr(kidof8,drhs,ida,jda,da,x,Ax,NAp,NAi,index1,limax8,0)
      print*, 'out of Umfpack'
C---------------------------------------------------
c      irefine=0       
c      enddo
C---------------------------------------------------
c     output the result
      open(2,file='../dat/result/cipsolution.dat')
      write(2,*)kdof
      do k=1,kbdof
        u(k)=-ub(k)
        write(2,*)-ub(k)
      enddo
      do k=1,kidof
        u(kbdof+k)=-x(k)
        write(2,*)-x(k)
      enddo
      close(2)
C---------------------------------------------------
      open(3,file='../dat/result/ltgmap.dat')
      write(3,*)nel
      do k1=1,nel
      do k2=1,4
        write(3,*)lgmap(k1,k2)
      enddo
      enddo
      close(3)
C---------------------------------------------------     
      open(1,file='../dat/refineinfo.dat')
      write(1,*)mx,my
      do j=1,my
       write(1,*)mxmax(j)
      enddo
      do mj=1,my
       write(1,*)(iflag(mi,mj),mi=1,mxmax(mj))
      enddo
      close(1)
      print*,'Wait computing the error'
C---------------------------------------------------      
C     计算误差
C     重写数值解格式
      call get_error(lgmap,nel,mxmax,i_fnum,mxy,iflag,
     &            u,kdof,mx,nz,p,symm,sigma1,sigma2,ineps)
C---------------------------------------------------

C---------------------------------------------------     
      stop
      end
