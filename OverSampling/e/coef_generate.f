	program coef_generate
C	to generate the permeability 
        implicit double precision(a-h,o-z)
	dimension perm(-511:2560,-511:2560)
	common/epsilon/eps_1,eps_2

c	nx=1024
c	ny=1024
	open(1,file='input.dat')
        read(1,*) levels,nxL,nxG,levelsglobal,ineps,mi,mj 
 
        read(1,*)
        read(1,*) xleft,yleft,xright,yright
	close(1)
	
 	eps_1=1.d0/ineps
	eps_2=eps_1
	nx=nxL*nxG
	ny=nx
	hx=(xright-xleft)/nx
	hy=(yright-yleft)/ny

        print*, 'eps_1=',eps_1, ' eps_2= ',eps_2
        if (nx .eq. 2048) then
                nxstart=-511
                nystart=-511
        else if (nx .eq. 1024) then
                nxstart=-255
                nystart=-255
        else if (nx .eq. 512 ) then
                 nxstart=-127
                 nystart=-127
        else if (nx .eq. 256 ) then
                 nxstart=-63
                 nystart=-63
        else if (nx .eq. 128) then
                 nxstart=-31
                 nystart=-31
        else if (nx .eq. 64) then
                 nxstart=-15
                 nystart=-15
        endif
 

       do i=nxstart,nx-nxstart+1
          do j=nystart,ny-nystart+1
 
             perm(i,j)=(coef(xleft+(i-1)*hx,yleft+(j-1)*hy)
     +		       +coef(xleft+(i-1)*hx,yleft+j*hy)
     +                 +coef(xleft+i*hx,yleft+(j-1)*hy)
     +		       +coef(xleft+i*hx,yleft+j*hy))/4.0d0
          end do
       end do
	
	nx1=nx-2*(nxstart-1)
	ny1=ny-2*(nystart-1)
	open(1,file='../dat/perms.dat')
	   write(1,*) nx1,ny1,nx1*ny1
	   do i=nxstart,nx-nxstart+1
	      write(1,*) (perm(i,j),j=nystart,ny-nystart+1 )
	   enddo
	close(1)

	end
