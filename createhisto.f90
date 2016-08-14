!--------------------------
program createhisto
!--------------------------
!
!   
!     Program to creat the data needed for histogram plot
!     (c) Sheng. Jin
!

IMPLICIT NONE

! Parameter readin
Logical                                  ::logbin       ! logarithm bin or not
INTEGER                                  ::nbin         ! number of bins
DOUBLE PRECISION                         ::xmin,xmax    ! The entire range of the bin
CHARACTER(LEN=200)                       ::filename


!parameters used to read the data file
INTEGER                                  ::i,itot
INTEGER                                  ::ii
INTEGER,PARAMETER                        ::MAXLINE=100000    
DOUBLE PRECISION                         ::bog       ! bog
DOUBLE PRECISION,DIMENSION(MAXLINE)      ::x_array   ! The VALUEs
INTEGER,DIMENSION(MAXLINE)               ::n_array   ! which bin each VALUE should be


! set the number of bins
INTEGER,dimension(:),allocatable                  ::xNnbin       ! Num 
DOUBLE PRECISION,dimension(:),allocatable         ::xNNnbin      ! Num/Ntot
DOUBLE PRECISION,dimension(:),allocatable         ::xmid         ! Mid-value
INTEGER                                           ::j


DOUBLE PRECISION,dimension(:),allocatable         ::x_grid       ! (Mid-value, Num/Ntot) , nbin
DOUBLE PRECISION                                  ::h_grid       ! grid size
INTEGER                                           ::k


!Logical,PARAMETER ::debug=.FALSE.
Logical,PARAMETER ::debug=.TRUE.


!read in initial conditions
OPEN(unit=3,file='createhisto.in')
READ(3,*)nbin
READ(3,*)logbin
READ(3,*)xmin
READ(3,*)xmax
READ(3,*)filename
!filename='RjReRj23Re23.dat'


IF(DEBUG)THEN
   WRITE(*,*)
   WRITE(*,*)'Parameters readin'
   WRITE(*,*)'nbin,logbin',nbin,logbin
   WRITE(*,*)'xmin,xmax  ',xmin,xmax
   WRITE(*,*)'filename   ',filename
   WRITE(*,*)
END IF


if(ALLOCATED(xNNnbin))then
     deallocate(xNNnbin)
end if
if(ALLOCATED(xNnbin))then
     deallocate(xNnbin)
end if
if(ALLOCATED(xmid))then
     deallocate(xmid)
end if
if(ALLOCATED(x_grid))then
     deallocate(x_grid)
end if

allocate(xNNnbin(nbin))
allocate(xNnbin(nbin))
allocate(xmid(nbin))
allocate(x_grid(nbin+1))



!set the bin grid
if(logbin)then
    h_grid=(log(xmax)-log(xmin))/real(nbin)  
else
    h_grid=((xmax-xmin)/nbin)
end if

IF(DEBUG)THEN
    write(*,*)'xmin,xmax,nbin',xmin,xmax,nbin
    write(*,*)'h_grid,logbin',h_grid,logbin
END IF

k=1
Do k=1,nbin+1
    if(logbin)then
        x_grid(k)=exp(log(xmin)+real(k-1)*h_grid) 
    else
        x_grid(k)=xmin+(k-1)*h_grid
    end if
    IF(DEBUG)THEN
         write(11,*),x_grid(k)
    END IF
end do

!stop

!read the data file and do statistics
open(unit=20,file=trim(filename),form='formatted',access='sequential',err=2001)
ii=0
DO i=1,MAXLINE
    !read
    read(20,fmt=*,end=2020),bog,bog,bog,x_array(i) 
    
    !find the bin where x_array(i) should be located
    CALL huntSM(x_grid,x_array(i),n_array(i),nbin+1)
    IF(DEBUG)THEN
        if ( (n_array(i).ne.0) .and. (n_array(i).ne.(nbin+1)) )then
            write(12,*)n_array(i),x_array(i),x_grid(n_array(i)),x_grid(n_array(i)+1)
        end if
    END IF

    ii=ii+1
END DO
2020 close(unit=20)
itot=ii
 

!do statistic and output
xNnbin=0
xNNnbin=0
DO j=1,itot
    if ( (n_array(j).ne.0) .and. (n_array(j).ne.(nbin+1)) )then
        xNnbin(n_array(j))=xNnbin(n_array(j))+1
    end if
END DO

DO i=1,nbin
    xNNnbin(i)=real(xNnbin(i))/real(itot)
    if(logbin)then
        xmid(i)=exp(log(x_grid(i))+h_grid/2.0)
!        x_grid(k)=exp(log(xmin)+real(k-1)*h_grid) 
    else
        xmid(i)=(x_grid(i)+x_grid(i+1))/2.0
    end if
END DO

DO i=1,nbin
    write(1,*),xmid(i), xNNnbin(i), xNnbin(i)
END DO


IF(DEBUG)THEN
     write(*,*)'itot',itot
     do i=1,itot
         write(14,*),x_array(i),n_array(i)
     end do
END IF


return
 
2001 write(*,*)'Error: file not find:',trim(filename)


END PROGRAM createhisto





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE huntSM(xx,x,jlo,N) 

IMPLICIT NONE 
INTEGER, INTENT(IN)                            :: N
INTEGER, INTENT(INOUT)                         :: jlo 
DOUBLE PRECISION, INTENT(IN)                   :: x 
DOUBLE PRECISION, DIMENSION(N), INTENT(IN)     :: xx 
!Given an array xx(1:N),and given a value x returns a va ue jlo such that x is between 
!xx(jlo) and xx(jlo+1).xx must be monotonic,either increasing or decreasing.
!jlo =0 or jlo = N is returned to indicate that x is out of range.jlo on input is taken as 
!the initial guess for jlo on output. 
INTEGER                                        :: inc,jhi,jm 
LOGICAL                                        :: ascnd 

 
ascnd = (xx(n) >= xx(1))                                 !True if ascending order of tab e,fa se otherwise. 
IF (jlo <= 0 .OR. jlo > n) THEN                          !Input guess not useful.Go immediately to bisection.
 jlo=0 
 jhi=n+1 
ELSE 
 inc=1                                                    !Set the hunting increment. 
 IF (x >= xx(jlo) .EQV. ascnd) THEN                       !Hunt up: 
   DO 
     jhi=jlo+inc 
     IF (jhi > n) THEN                                    ! Done hunting,since off end of table 
       jhi=n+1 
       EXIT 
     ELSE 
       IF (x < xx(jhi) .EQV. ascnd) EXIT 
       jlo=jhi                                           !Not done hunting, 
       inc=inc+inc                                       !so double the increment 
    END IF 
  END DO                                                 !and try again. 
ELSE                                                     !Hunt down: 
   jhi=jlo 
   DO 
     jlo=jhi-inc 
     IF (jlo < 1) THEN                                   !Done hunting,since off end of table
          jlo=0 
          EXIT 
     ELSE 
       IF (x >= xx(jlo) .EQV. ascnd) EXIT 
       jhi=jlo                                           !Not done hunting, 
       inc=inc+inc                                       !so double the increment 
     END IF 
  END DO                                                 ! and try again. 
 END IF 
END IF                                                   !Done hunting,va ue bracketed.
DO                                                       !Hunt is done,so begin the finnal bisection phase: 
 IF (jhi-jlo <= 1) THEN 
   IF (x == xx(n)) jlo=n-1 
   IF (x == xx(1)) jlo=1 
   EXIT 
 ELSE 
   jm=(jhi+jlo)/2 
   IF (x >= xx(jm) .EQV. ascnd) THEN 
     jlo=jm 
   ELSE 
     jhi=jm 
  END IF 
 END IF 
END DO 
END SUBROUTINE huntSM

