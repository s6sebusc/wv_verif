subroutine vgram_f( nx, ny, y ,r ,newr, newl)
    implicit none
    integer :: nx, ny, i, j, r
    real(kind=8), dimension((r+1),(r+1)) :: newr, newl
    real(kind=8), dimension(nx,ny) :: y    
    newr(:,:)=0
    newl(:,:)=0
      do i = 0,r
        do j = 0,r
        if(j**2+i**2 <= r**2) then 
          newr((i+1),(j+1))=0.5*sum((abs(y((i + 1):nx,(j + 1):ny) - y(1:(nx-i),1:(ny-j)))) )&
          /((nx-i)*(ny-j))
          newl((i+1),(j+1))=0.5*sum((abs(y(1:(nx-i),(j + 1):ny) - y((i + 1):nx,1:(ny-j)))) )&
          /((nx-i)*(ny-j))
          end if
       end do
     end do
end subroutine vgram_f


