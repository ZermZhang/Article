      subroutine inverse(A,n)
      implicit double precision(a-h,o-z) 
      dimension A(n,n),IP(n) 
C      dimension A(3,3),IP(3)
      do 100 k=1,n !step1
	ir=k    !step2
	pmax=dabs(A(ir,k))
	do i=k+1,n
          if(pmax.lt.dabs(A(i,k))) then
            pmax=dabs(A(i,k))
            ir=i
          endif
	end do
        IP(k)=ir
        if(A(IP(k),k).eq.0.d0) then  !step3
          print*,"A is singular!"
          stop
        endif
        if(IP(k).ne.k) then!step4
          do i=1,n
            u=A(IP(k),i)
	    A(IP(k),i)=A(k,i)
	    A(k,i)=u
	  end do
        endif
        A(k,k)=1.0d0/A(k,k) !step5
        do i=1,n !step6
          if(i.ne.k) then
	    A(i,k)=-A(k,k)*A(i,k)
	  endif
        enddo
        do i=1,n 	!step7
	  if(i.ne.k) then
            do j=1,n
	      if(j.ne.k) then
                A(i,j)=A(i,j)+A(i,k)*A(k,j)
              endif
            end do
          endif
        end do
        do j=1,n
	  if(j.ne.k) then
	    A(k,j)=A(k,k)*A(k,j)
          endif
        enddo
100   continue           

      do k=n,1,-1	!step9
	if(k.ne.IP(k)) then
	  do i=1,n
	    u=A(i,IP(k))
	    A(i,IP(k))=A(i,k)
	    A(i,k)=u
	  enddo
        endif
      enddo
      return
      end
