! Compute the time informational entropy on fully connected networks
! The N1 and N0 frozen nodes are connected to all lattice points at the beginning
! N1 and N0 connections disappear biased by neighbour influence
!
! Marcus A.M. de Aguiar and Gabriella D. Franco - Modified - 18/01/2022

program voter_dyn_ext

IMPLICIT REAL*8 (A-H,O-Z)
REAL*8 aux,aux2,anr,annb,prba,auxi,annbi, d2
REAL*8 ai, tmax, p, asym, lambda, np, rn, r2, n0f, n1f
INTEGER, ALLOCATABLE :: x(:),b(:,:),pr(:),viz(:), bi(:,:), a(:,:), nn(:)
REAL*8, ALLOCATABLE :: n0(:),n1(:)
INTEGER  i,j,k,h,u,ik,imk,tot,kn,knviz,l,ii,deg, clusters, clustersize, counts, counts_nets, tau
INTEGER  n,nt,nm,nnb, ind, samples, nsamples, del_m
INTEGER  inew, jnew, ri, pp, q, column, line, seq, rew, icheck, newneigh, i1
INTEGER :: iseed(12)
CHARACTER*70 filename,filename2,filename3, filename4
CHARACTER*3 strings,strdeg, strclu, strsample


OPEN(UNIT=50,FILE='seed.in',STATUS='OLD')
READ(50,*) iseed
CLOSE(50)
CALL RANDOM_SEED(put=iseed)


OPEN(UNIT=8,FILE='input.in',STATUS='OLD',POSITION='REWIND')
	!model parameters
	READ(8,*)n            !population size
	READ(8,*)deg          !network degree = 2*deg
	READ(8,*)np           !total value of external influence
	READ(8,*)tau          !equilibrium time
	READ(8,*)nm           !number of measures to acquire (lines in the output file)
	READ(8,*)del_m        !interval between measures
	READ(8,*)lambda       !step used in external influence changes
CLOSE(8)


an = dfloat(n)
anb = dfloat(nm)

nmax = tau+del_m*nm

nsamples = 1 !if more samples are wanted



ALLOCATE (x(n),pr(0:n),b(n,n-1),n0(n),n1(n),viz(n),bi(n,n),a(n,n),nn(n))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generating topology
! Other topologies can be used. See commented sections at the end


!Generate Ring Network with degree 2*deg
b=0
do i=1,n
    do j=1,deg
        ik = i + j
        imk = i - j
        if(ik > n) ik = ik - n
        if(imk < 1) imk = imk + n
        b(i,j) = imk
        b(i,j+deg) = ik
    end do
end do

!Output Files
filename = 'ext_corr.dat'
OPEN(UNIT=58,FILE=filename,STATUS='UNKNOWN')

filename2 = 'network.dat'
OPEN(UNIT=20,FILE=filename2,STATUS='UNKNOWN')

filename3 = 'ext_sum.dat'
OPEN(UNIT=37,FILE=filename3,STATUS='UNKNOWN')

! INITIAL CONFIGURATION
pr = 0
prba = 0.0d0
! dynamical probability
p = 0


cluster=0

simulations: DO samples = 1,nsamples
    x = 0 ! Initial configuration


    ! Initializing network state randomly :: Marcus
    initial_state: DO i=1,n
        call random_number(aux)
        if(aux > 0.5) then
            x(i) = 1
            !----Initializing external vector
            !----n1(j):=number of external
            !---influencers acting on a node j
            call random_number(aux2)
            n1(i) = aux2*np !initial bias :: Gabi
            n0(i) = np - n1(i)
        else
            call random_number(aux2)
            n1(i) = aux2*np
            n0(i) = np - n1(i)
        end if
    END DO initial_state
    !write(37,*) x
!    pause

    nf = 1 !selecting when to acquire measures
    counts = 0
    counts_nets = 0

    ! DYNAMICS
    dynamics: DO i=1,nmax
        if(i > tau) nf = mod(i,del_m)
        call random_number(aux)
        j = int(aux*n) + 1 ! selecting focal node
        call random_number(aux)

        kn = count(b(j,:)>0) !  counts the number of neighbors j has inside the dynamic network
        annb = dfloat(kn)
        if(aux > p) then
            call random_number(aux)
            aux2 = aux*(annb+n0(j)+n1(j))
            if(aux2 > annb + n0(j)) then
                x(j) = 1
            else if(aux2 > annb) then
                x(j) = 0
            else
                call random_number(aux) ! this new random variable was  missing
                aux2 = aux*(annb)
                k = int(aux2) + 1
                x(j) = x(b(j,k))
            end if
        end if


        asym = 0
        do h=1,kn
            asym = asym + x(b(j,h))/(annb)
        end do
        !print *, asym

        n0f = n0(j)
        n1f = n1(j)
        if(asym > 0.5d0 .AND. n0(j)>=lambda) then
            n0(j) = n0(j) - lambda
            n1(j) = n1(j) + lambda
        else if(asym < 0.5d0 .AND. n1(j)>=lambda) then
            n1(j) = n1(j) - lambda
            n0(j) = n0(j) + lambda
        else if(asym==0.5d0 .AND. n0(j)>=lambda .AND. n1(j)>=lambda) then
            CALL random_number(aux)
            if(int(2*aux)==1)then
                n1(j) = n1(j) + lambda
                n0(j) = n0(j) - lambda
            else
                n1(j) = n1(j) - lambda
                n0(j) = n0(j) + lambda
            end if
        end if

        if(counts /= nm .AND. nf==0) then
            write(58,*) DOT_PRODUCT(n0,(1-x))/sum(n0)
            write(58,*) DOT_PRODUCT(n1,x)/sum(n1)
            write(37,*) sum(n0), sum(n1)
            counts = counts + 1
            clusters=0
            clustersize = 1
            DO ii=2,n
                if(x(ii) == x(ii-1)) then
                    clustersize = clustersize+1
                else
                    write(20,*) clustersize
                    clustersize = 1

                    clusters = clusters+1
                    write(20,*) x(ii-1)
                end if
            END DO
            write(20,*) clustersize
            clusters=clusters+1
            write(20,*) x(n)
       end if

        if(nf == 0) then
            tot = sum(x)  ! count number of states 1
            pr(tot) = pr(tot) + 1
        end if
    END DO dynamics

        DO ii=0,n
            prba = dfloat(pr(ii))/(dfloat(nsamples)*anb)
        END DO


        print *, 'simulation ok'
        print *, sum(n0), sum(n1)

END DO simulations

CLOSE(58)

CLOSE(37)

CLOSE(20)

CALL RANDOM_SEED(get=iseed)
OPEN(UNIT=50,FILE='seed.in',STATUS='OLD', POSITION='REWIND')
WRITE (50,*) iseed
close(50)

end program voter_dyn_ext


FUNCTION seq(i,j,l) result(k)
  INTEGER, intent(in) :: i,j
  INTEGER :: k
  k = l*(i-1)+j
END FUNCTION seq


FUNCTION line(k,l) result(i)
  INTEGER, intent(in) :: k
  INTEGER :: i
  i = INT((k-1)/l)+1
END FUNCTION line

FUNCTION column(k,l) result(j)
  INTEGER, intent(in) :: k
  INTEGER :: j
  j = k-l*INT((k-1)/l)
END FUNCTION column

SUBROUTINE FC_gen(n,b)
    INTEGER, intent(in) :: n
    INTEGER :: i,j
    INTEGER :: b(n,n-1)
    b=0
    do i=1,n
        do j=i+1,n
            b(i,j-1) = j
            b(j,i) = i
        end do
    end do

END
! Generate Totally Connected Network
! b(i,j) = j-th neighbor of node i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE NUMBSTR(ID,NUMBER,STR)
CHARACTER*(*) STR
INTEGER*4 ID,NUMBER
CHARACTER*1 B
INTEGER*4 IA0,NN,II,IT
IA0 = ICHAR('0')
NN = NUMBER
DO II=1,ID
J = ID + 1 - II
IT = MOD(NN,10)
B = CHAR(IA0 + IT)
STR(J:J) = B
NN = (NN - IT)/10
END DO
RETURN
END
