module parameters
  IMPLICIT NONE
  !-- Parameters -------------------------------------------------
  integer, parameter :: D=2           !D=2 or D=3  
  integer, parameter :: ScaleNum=64    !number of scales
  integer, parameter :: ExtMomNum=32
  integer, parameter :: AngleNum=32
  integer, parameter :: kNum=512       !k bins of Green's function
  integer, parameter            :: MaxOrder=4           ! Max diagram order

  double precision   :: UVScale=100.0     !the upper bound of energy scale
  double precision   :: MaxExtMom=4.0     !the upper bound of energy scale
  double precision   :: DeltaW
  double precision   :: DeltaQ

  integer            :: PID           ! the ID of this job
  integer            :: Order
  double precision   :: Mass2         ! square mass
  double precision   :: Beta, Mu, Rs, Kf, Ef
  double precision   :: BareCoupling         ! bare coupling
  integer            :: Seed          ! random-number seed

  !-- Markov Chain ----------------------------------------------
  double precision                        :: Step        ! a counter to keep track of the current step number
  integer, parameter                      :: UpdateNum=5 ! number of updates
  double precision, dimension(UpdateNum)  :: PropStep, AcceptStep
  double precision, dimension(0:MaxOrder) :: ReWeightFactor       !reweightfactor for each order

  integer                                 :: CurrOrder, CurrScale, CurrIRScale
  integer                                 :: CurrExtMom  !external momentum for self energy
  double complex                          :: CurrWeight
  double precision, dimension(D+1, MaxOrder+3)  :: LoopMom ! values to attach to each loop basis
  double precision, dimension(D+1, ExtMomNum)        :: ExtMomMesh ! external momentum
  double precision, dimension(D+1)          :: Mom0 ! values to attach to each loop basis

  !--- Measurement ------------------------------------------------
  double precision, dimension(ScaleNum)       :: ScaleTable, dScaleTable
  double precision, dimension(MaxOrder, ScaleNum, AngleNum, ExtMomNum)       :: DiffVer
  double precision       :: Norm

  !-- Diagram Elements  --------------------------------------------
  double precision, dimension(ScaleNum, AngleNum, ExtMomNum)       :: EffVer

  !-- common parameters and variables ------------------------------
  ! THIS IS PROJECT-INDEPENDENT 
  integer, parameter          :: UP=1, DOWN=0
  double precision, parameter :: tm32   = 1.d0/(2.d0**32.d0)
  double precision, parameter :: eps    = 1.d-14            ! very small number
  integer,          parameter :: Mxint  = 2147483647        ! maximum integer
  integer,          parameter :: Mnint  =-2147483647        ! minimum integer
  double precision, parameter :: pi=3.1415926
end module

INCLUDE "rng.f90"

program main
    use mt19937
    use parameters
    implicit none
    integer :: iBlck, iInner, AnnealStep, TotalStep
    double precision :: x
  
    print *, 'Beta, rs, mass2, Order, TotalStep(*1e6), Seed, PID'
    read(*,*) Beta, Rs, Mass2, Order, TotalStep, Seed, PID

    ! For a given order, the bigger factor, the more accurate result 

    call sgrnd(Seed) 
  
    print *, "Initializing ..."
    call Initialize()
    print *, "Initialized!"
  
    call Test() !call test first to make sure all subroutines work as expected
  
    AnnealStep=4
    Step=0.0
    print *, "Start simulation..."
    do iBlck=1,TotalStep
      do iInner=1,1000000
        Step=Step+1.0
        x=grnd()
        if (x<1.0/UpdateNum) then
          call ChangeOrder()
        else if (x<2.0/UpdateNum) then
          call ChangeMom()
        else if (x<3.0/UpdateNum) then
          call ChangeExtQ()
        else if (x<4.0/UpdateNum) then
          call ChangeExtK()
        ! else if (x<2.0/UpdateNum) then
        !   call ChangeScale()
        endif
        ! if(abs(LoopMom(D+1, 4))>UVScale) then
        !   print *, x, Step
        !   stop
        ! endif
        !if(mod(int(Step), 4)==0) call Measure()
        call Measure()

      enddo

      ! call DynamicTest()
      ! call SolveBetaFunc()

      ! if (iBlck==AnnealStep) then
      !   AnnealStep=AnnealStep*2
      !   CurrIRScale=CurrIRScale/2
      ! endif

      if(mod(iBlck, 10)==0) then
        !!!!!!!!!!!!  Print Info and Save Data !!!!!!!!!!!!!!!!!!!!!!!
        print *, "freq:", LoopMom(D+1, 4)
        print *, "mom:", norm2(LoopMom(1:D, 4))
        print *, "Weight: ", CurrWeight
        ! print *, "Norm:", Norm
        ! print *, DiffVer(1, 1, 1, :)
        write(*,*) 
        print *, iBlck, "million steps"
        write(*,"(A20)") "Accept Ratio: "
        write(*,"(A16, f8.3)") "Increase Order:", AcceptStep(1)/PropStep(1)
        write(*,"(A16, f8.3)") "Decrease Order:", AcceptStep(2)/PropStep(2)
        write(*,"(A16, f8.3)") "Change Mom:", AcceptStep(3)/PropStep(3)
        write(*,"(A16, f8.3)") "Change ExtQ:", AcceptStep(4)/PropStep(4)
        write(*,"(A16, f8.3)") "Change ExtK:", AcceptStep(5)/PropStep(5)
        ! write(*,"(A16, f8.3)") "Change Scale:", AcceptStep(5)/PropStep(5)
        ! write(*, *) "coupling: ", DiffVer(CurrScale)/Norm
        ! write(*, *) "coupling: ", EffVer
        ! write(*, *) "coupling: ", DiffVer/Norm
        ! write(*, *) "mu: ", DiffMu/Norm
        ! write(*, *) "mu: ", EffMu
      endif

      if (mod(iBlck, 100)==10)  then
        call SaveToDisk()
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end do

    call SaveToDisk()
    print *, "End simulation."
  
    CONTAINS

    subroutine Initialize()
      implicit none
      integer :: i, j
      double precision :: kamp
  ! For a given order, the bigger factor, the more accurate result 
      Kf=sqrt(2.0)/Rs
      Ef=Kf*Kf
      Mu=Ef
      Beta=Beta/Ef
      MaxExtMom=MaxExtMom*Kf
      UVScale=UVScale*Kf
      DeltaW=2.0*pi/Beta

      ReWeightFactor(0:2)=(/100.0,1.0,20.0/)

      PropStep=0.0
      AcceptStep=0.0

      do i=1, ScaleNum
        ScaleTable(i)=(i-1)*1.0/ScaleNum*UVScale
      end do
      ScaleTable(1)=1.0e-6
      do i=1, ScaleNum-1
        dScaleTable(i)=ScaleTable(i+1)-ScaleTable(i)
      end do

      ExtMomMesh=0.0
      DeltaQ=MaxExtMom/ExtMomNum
      do j=1, ExtMomNum
        !1 <--> DeltaK/2, kNum <--> kMax-DeltaK/2
        kamp=(j-0.5)*DeltaQ
        ExtMomMesh(1, j)=kamp
      enddo
      Mom0=0.0
      Mom0(1)=Kf
      Mom0(D+1)=DeltaW*0.5

      EffVer=1.0

      DiffVer=0.0
      Norm=1.0e-10

      CurrScale=4
      CurrIRScale=ScaleNum/2
      CurrOrder=0
      CurrExtMom=1

      LoopMom=0.0
      do i=4, 3+Order
        LoopMom(1, i)=Kf
        LoopMom(D+1, i)=DeltaW*0.5
      enddo
      LoopMom(:, 1)=ExtMomMesh(:, CurrExtMom)
      LoopMom(:, 2)=Mom0
      LoopMom(:, 3)=Mom0

      CurrWeight=CalcWeight(CurrOrder)
    end subroutine

    subroutine Test()
      implicit none
      return
    end subroutine

    subroutine DynamicTest()
      implicit none
      return
    end subroutine

    subroutine Measure()
      implicit none
      double precision :: Factor

      ! if(CurrIRScale>=CurrScale) return

      Factor=real(CurrWeight)/abs(CurrWeight)/ReWeightFactor(CurrOrder)
      if(CurrOrder==0) then
          Norm=Norm+Factor
      else
          DiffVer(CurrOrder, 1, 1, CurrExtMom)=DiffVer(CurrOrder, 1, 1, CurrExtMom)+Factor
      endif
      return
    end subroutine
    
    subroutine SaveToDisk()
      implicit none
      integer :: i, ref, j, o
      double precision :: Obs
      character*10 :: ID
      character*10 :: order_str
      character*20 :: filename
      !Save Polarization to disk

      ! do o=1, Order
        write(ID, '(i10)') PID
        write(order_str, '(i10)') o
        filename="Data/polar"//trim(adjustl(order_str))//"_pid"//trim(adjustl(ID))//".dat"
        write(*,*) "Save to disk ..."
        open(100, status="replace", file=trim(filename))
        write(100, *) "#", Step
        do i=1, ExtMomNum
            Obs = DiffVer(1, 1, 1, i)/Norm
            write(100, *) norm2(ExtMomMesh(1:D, i)), Obs
        enddo
        close(100)
      ! enddo

      return
    end subroutine

    subroutine SolveBetaFunc()
      implicit none
      integer :: i, start, end
      double precision :: dg, dMu
    end subroutine

    double complex function CalcWeight(NewOrder)
      !calculate the weight for ALL diagrams in a given sector
      implicit none
      integer :: NewOrder
      if(NewOrder==0) then 
        ! CalcWeight=1.0/(ExtMomMesh(1, CurrExtMom)**2+1.0)
        CalcWeight=1.0
      else if(NewOrder==1) then
        CalcWeight=Ver4Loop1(NewOrder, Mom0, Mom0, ExtMomMesh(:, CurrExtMom))
      endif
      ! print *, CalcWeight
      CalcWeight=CalcWeight*ReWeightFactor(NewOrder)
      return
    end function CalcWeight
    
    subroutine ChangeOrder()
      !increase diagram order by one/change normalization diagram to physical diagram
      implicit none
      double precision :: R, Prop
      double complex :: Weight
      integer :: NewOrder, Index
      if(grnd()<0.5) then
        ! increase order
        if (CurrOrder==Order) return
        Index=1
        NewOrder=CurrOrder+1
        call CreateMom(LoopMom(:, CurrOrder+4), Prop)
        if(Prop<=0.0) return
      else
        if (CurrOrder==0) return
        Index=2
        NewOrder=CurrOrder-1
        call RemoveMom(LoopMom(:, CurrOrder+3), Prop)
        if(Prop<=0.0) return
      endif

      PropStep(Index)=PropStep(Index)+1.0
      Weight = CalcWeight(NewOrder)
      ! print *, Weight, NewOrder, Index
      R=abs(Weight)/abs(CurrWeight)*Prop
      if(grnd()<R) then
        AcceptStep(Index)=AcceptStep(Index)+1.0
        CurrWeight=Weight
        CurrOrder=NewOrder
      endif
      return
    end subroutine

    subroutine ChangeMom()
      !randomly choose a vertex, change the space variable
      implicit none
      double precision, dimension(D+1) :: NewMom, OldMom
      double precision :: prop, R
      double complex :: Weight
      integer :: Num
  
      if(Order==0) return
      Num = int( CurrOrder*grnd() ) + 4
  
      PropStep(3) = PropStep(3) + 1.0
      OldMom = LoopMom(:, Num)
  
      call ShiftMom(OldMom, NewMom, prop)
  
      if(prop<=0.0) return
      LoopMom(:,Num)=NewMom
  
      Weight = CalcWeight(CurrOrder)
      R = prop*abs(Weight)/abs(CurrWeight)
  
      if(grnd()<R) then
        AcceptStep(3) = AcceptStep(3)+1.0
        CurrWeight = Weight
      else
        LoopMom(:,Num)=OldMom
      endif
  
      return
    end subroutine

    subroutine ChangeExtMom()
      !randomly choose a vertex, change the space variable
      implicit none
      double precision :: prop, R
      double complex :: Weight
      integer :: OldExtMom
      integer, parameter :: Delta=3

      OldExtMom=CurrExtMom
      if(grnd()<0.5) then
        CurrExtMom=CurrExtMom+int(grnd()*Delta)
      else
        CurrExtMom=CurrExtMom-int(grnd()*Delta)
      endif
      if(CurrExtMom>ExtMomNum .or. CurrExtMom<1) then
        CurrExtMom=OldExtMom
        return
      endif

      LoopMom(:, 1)=ExtMomMesh(:, CurrExtMom)
  
      PropStep(4) = PropStep(4) + 1.0
      Weight = CalcWeight(CurrOrder)
      R = abs(Weight)/abs(CurrWeight)
  
      if(grnd()<R) then
        AcceptStep(4) = AcceptStep(4)+1.0
        CurrWeight = Weight
        ! call UpdateState()
      else
        CurrExtMom=OldExtMom
        LoopMom(:, 1)=ExtMomMesh(:, CurrExtMom)
      endif
      return
    end subroutine

    subroutine ChangeExtK()
      !randomly choose a vertex, change the space variable
      implicit none
      double precision :: R, theta
      double complex :: Weight
      double precision, dimension(D+1) :: OldK
      integer :: Num
      integer, parameter :: Delta=3

      theta=grnd()*2.0*PI
      if(grnd()<0.5) then
        Num=2
      else
        Num=3
      endif

      OldK=LoopMom(:, Num)

      LoopMom(1, Num)=Kf*cos(theta)
      LoopMom(2, Num)=Kf*sin(theta)
  
      PropStep(5) = PropStep(5) + 1.0
      Weight = CalcWeight(CurrOrder)
      R = abs(Weight)/abs(CurrWeight)
  
      if(grnd()<R) then
        AcceptStep(5) = AcceptStep(5)+1.0
        CurrWeight = Weight
        ! call UpdateState()
      else
        LoopMom(:, Num)=OldK
      endif
      return
    end subroutine

    subroutine ChangeScale()
      implicit none
      !TODO: don't forget to change all LoopMom with scale!
      integer :: OldScale
      double complex :: Weight
      double precision :: R

      OldScale=CurrScale
      if(grnd()<0.5) then
        CurrScale=CurrScale-1
      else
        CurrScale=CurrScale+1
      endif

      if(CurrScale<1 .or. CurrScale>ScaleNum) then
        CurrScale=OldScale
        return
      endif

      PropStep(3) = PropStep(3) + 1.0

      Weight = CalcWeight(CurrOrder)
      R = abs(Weight)/abs(CurrWeight)
      if(grnd()<R) then
        AcceptStep(3) = AcceptStep(3)+1.0
        CurrWeight = Weight
      else
        CurrScale=OldScale
      endif

      return
    end subroutine

    subroutine CreateMom(New, Prop)
      implicit none 
      double precision, dimension(D+1) :: New
      double precision :: Prop, dK, Kamp, theta, phi
      ! dK=1.0*CurrScale
      dK=1.0*ScaleTable(CurrScale)
      Kamp=Kf+(grnd()-0.5)*2.0*dK
      if(Kamp<=0.0) then
        Prop=-1.0
      endif

      phi=2.0*pi*grnd()
      New(1)=kamp*cos(phi)
      New(2)=kamp*sin(phi)
      Prop=2.0*dK*2.0*pi*kamp

      !!!! the Hard way  !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Prop=-1.0
      ! if(Kamp<=0.0) return
      ! phi=2.0*pi*grnd()
      ! theta=pi*grnd()
      ! if(theta==0.0) return
      ! New(1)=Kamp*sin(theta)*cos(phi)
      ! New(2)=Kamp*sin(theta)*sin(phi)
      ! New(3)=Kamp*cos(theta)
      ! Prop=2.0*dK*2.0*pi*pi*sin(theta)*Kamp**(D-1)

      !!! Simple way  !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do i=1, D
      !   NewMom(i)=kF*(grnd()-0.5)*2
      ! enddo
      ! Prop=Beta*(2.0*kF)**D
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! frequency
      New(D+1)=(int((grnd()-0.5)*20)+0.5)*DeltaW
      Prop=Prop*20.0;

    end subroutine

    subroutine RemoveMom(Old, Prop)
      implicit none
      double precision, dimension(D+1) :: Old
      double precision :: Prop, dK, Kamp, SinTheta

      dK=1.0*ScaleTable(CurrScale)
      Kamp=norm2(Old(1:D))
      if(Kamp<Kf-dK .or. Kamp>Kf+dK) then
        Prop=-1.0
        return
      endif
      Prop=1.0/(2.0*dK*2.0*pi*Kamp)

      !Get proper K proposed probability
      !!!!!!! Hard way !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! SinTheta=norm2(Old(1:2))/Kamp
      ! if(SinTheta==0.0) then
      !   Prop=-1.0
      !   return
      ! endif
      ! Prop=1.0/(2.0*dK*2.0*pi*pi*SinTheta*Kamp**(D-1))

      !!!!!!! Simple way !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do i=1, D
      !   if(abs(LoopMom(i, LoopNum(CurrOrder)))>kF) return
      ! enddo
      ! Prop=1.0/(Beta*(2.0*kF)**D)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(Old(D+1)>10.5*DeltaW .or. Old(D+1)<-9.5*DeltaW) then
        Prop=-1.0
        return
      endif
      Prop=Prop/20.0

    end subroutine

    
    subroutine ShiftMom(old, new, prop)
      implicit none
      double precision, dimension(D+1) :: old, new
      integer :: Num
      double precision :: ratio, prop, k, k_new
      double precision :: lambda
      x=grnd()
      if(x<1.0/3.0) then
          new=old
          Num=int(grnd()*D)+1
          new(Num)=new(Num)+sign(grnd(), grnd()-0.5)
          prop=1.0
      else if(x<2.0/3.0) then
          k=norm2(old(1:D)) !sqrt(m_x^2+m_y^2+m_z^2)
          if(k==0.0)then
            prop=-1.0
            return
          endif
          lambda=1.5
          k_new=k/lambda+grnd()*(lambda-1.0/lambda)*k
          ratio=k_new/k
          new(1:D)=old(1:D)*ratio
          if(D==2) then
            prop=1.0 ! prop=k_old/k_new
          else
            prop=k_new/k ! prop=k_old/k_new
          endif
      else
          new(1:D)=-old(1:D)
          prop=1.0
      endif
      if(grnd()<0.5) then
        new(D+1)=old(D+1)+DeltaW
      else
        new(D+1)=old(D+1)-DeltaW
      endif
    end subroutine

    double complex function Green(Mom, Scale, g_type)
    !dimensionless green's function
      implicit none
      double precision :: kk, freq
      integer :: g_type, Scale
      double precision, dimension(D+1) :: Mom
      kk=norm2(Mom(1:D))
      freq=Mom(D+1)
      ! print *, "green", freq/DeltaW, freq, UVScale
      if(abs(freq)<UVScale) then
        Green=-1.0/(dcmplx(0.d0, freq)-kk*kk+Mu)
        ! Green=-1.0/(freq*freq+kk*kk+Mu)
      else
        Green=0.0
      endif
      ! Green=-1.0/(dcmplx(0.d0, freq)-kk*kk+Mu)
      ! print *, "green got", Green
      return
    end function

    double precision function Inter(InL, InR, Q)
      implicit none
      double precision, dimension(D+1) :: InL, InR, Q, ExtQ
      ExtQ=InL-InR-Q
      Inter=8.d0*PI*(1.0/(norm2(Q(1:D))**2+Mass2)-1.0/(norm2(ExtQ(1:D))**2+Mass2))
      return
    end function
    
    double complex function Ver4Loop1(Order, InL, InR, Q)
      implicit none
      integer :: Order
      double precision, dimension(D+1) :: InL, InR, Q, kL2R, kR2L 
      double precision :: k, freq
      kL2R=LoopMom(:, 4)
      kR2L=kL2R-Q
      Ver4Loop1=Green(kL2R, CurrScale, 0)*Green(kR2L, CurrScale, 0)
      ! Ver4Loop1=Ver4Loop1-Green(LoopMom(:, 4), CurrScale, 0)*Green(LoopMom(:, 4), CurrScale, 0)
      Ver4Loop1=Ver4Loop1*(Inter(InL, kR2L, Q)*Inter(kL2R, InR, Q)-(8.d0*PI/(norm2(Q(1:D))**2+Mass2))**2)
      ! Ver4Loop1=Ver4Loop1*(Inter(InL, kR2L, Q)*Inter(kL2R, InR, Q))

      ! k=norm2(LoopMom(1:D,4))
      ! freq=LoopMom(D+1, 4)

      ! Ver4Loop1=cmplx(exp(-k)*exp(-freq*freq/4.0)/(norm2(Q)**2+1), 0.0)

      ! Ver4Loop1=Ver4Loop1-Green(LoopMom(:, 4), CurrScale, 0)*Green(LoopMom(:, 4), CurrScale, 0)
      ! print *, LoopMom(D+1, 4), Ver4Loop1, CurrWeight
      return
    end function

end program main
