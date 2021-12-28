module ElasticFunctionsByIgor
implicit none  
contains 
subroutine MeshCreate()! creating X(i) Y(i) and W(i,j) for every i,j
 
end subroutine MeshCreate
double precision function w_An(x,y,R) !function w(x,y,t)
    double precision :: x, y, R, radius2, temp 
    radius2 = (x*x+y*y)
    temp =R*R-radius2
    w_An =sqrt( temp)
end function w_An

 double precision function f(x, y, x1, y1) !denominator of integral 
    double precision :: x, y, x1, y1,radius, temp
    radius= (x-x1)**2+(y-y1)**2
    temp=(1/radius)**3
    f=sqrt(temp)
    end function f

 double precision function  Gaus4points(x1,y1,x2,y2,h) !4 point Gauss method for denominator
! TODO: Is it possible to write Gauss in general case? not 4 points only. (Not inportant)

    double precision :: temp = 0
    double precision :: Root = 1/sqrt(3.0)
    double precision :: xi1, eta1, xi2, eta2, x1, y1, x2, y2,h
    xi1=-h*root/2 + x1
    eta1=-h*root/2 + y1
    xi2=h*root/2 + x1
    eta2=h*root/2 + y1
    temp = f( xi1, eta1,x2,y2)+f(xi1, eta2,x2,y2)+f( xi2, eta1,x2,y2)+f(xi2, eta2,x2,y2)
    Gaus4points = temp *h*h/4

 end function Gaus4points

 double precision function Get_p(k,l, R, h, n, X, Y,  W)! Evaluating of main integral as sum
    double precision :: R, h, radius, temp, singular, regular   
    integer :: n, i, j, k, l                   
    double precision, allocatable :: X(:), Y(:),  W(:,:) 
    temp = 0
    singular=0
    regular=0
    do i = 1, n
        do j = 1, n
            if ( i == k .and. j==l) then
                singular =-8*sqrt(2.0)*W(i,j)/h !singular part of the integral
                
                temp=temp+singular
               
            else
                regular= Gaus4points(X(i),Y(j),X(k),Y(l),h)* W(i,j)! regular part of the integral
                temp=temp + regular
            end if 
       end do
    end do
    Get_p = temp

 end function Get_p

subroutine Pressure(R, h, n, X, Y, P, W)! geting pressure and outputing in file
    double precision :: R, h, radius       
    integer :: n, i, j, k, l                   
    double precision, allocatable :: X(:), Y(:), P(:,:),  W(:,:) 
      do k = 1, n
        do l = 1, n
            P(k,l)=Get_p(k,l, R, h, n, X, Y, W)     
            
        enddo

    enddo
end subroutine Pressure

Subroutine MaxErr()     
   
end subroutine MaxErr

subroutine FileOutput(X,Y,P, n,h)
    double precision:: h
    double precision, allocatable :: P(:,:),X(:),Y(:)
    integer ::  k, l, n
        do k = 1, n
            do l = 1, n
                write (1,*) X(k)," ",Y(l)," ",P(k,l)       
            enddo
            write (1,*) "*"

    enddo
   
    
end subroutine FileOutput

 end module ElasticFunctionsByIgor

 program main
    use ElasticFunctionsByIgor
    double precision :: R, h, radius,temp       
    integer :: n, i, j                 
    double precision, allocatable :: X(:), Y(:), P(:,:),  W(:,:) ,P_EX(:,:)
    Read(*,*) R,n
    open(1, file ="out.txt")
    h=2*R/n
    allocate (X(n))
    allocate (Y(n))
    allocate (W(n,n))
    allocate (P(n,n))
    allocate (P_EX(n,n))

    do i = 1, n
        X(i)=-R+h/2+(i-1)*h
        Y(i)=-R+h/2+(i-1)*h
    enddo

    do i = 1, n
        do j= 1, n
        radius = sqrt((X(i)*X(i)+Y(j)*Y(j)))
     
            if (radius <=R) then
                P_EX(i,j)=1
                !W(i,j)=w_An(X(i),Y(j),R)
                W(i,j)=1
            else 
                P_EX(i,j)=0
                W(i,j) = 0
            endif
            
        enddo
     enddo
   !call MeshCreate     
    call Pressure(R, h, n, X, Y, P, W)
   
    call FileOutput(X,Y,P,n,h)
    !call MaxErr(P,P_EX)
    temp=maxval(abs(P_EX-P))
    write(1,*) "MaxErr= ", temp
    deallocate(X)
    deallocate(Y)
    deallocate(W)
    deallocate(P)
    deallocate(P_EX)
end program main