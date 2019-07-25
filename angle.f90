module angle2DUtility
    implicit none
    double precision, parameter, private :: pi=3.1415926536
    integer, parameter, private :: D=2
    double precision, parameter, private :: eps    = 1.d-14            ! very small number

    private :: TestAngle2D, TestAngleIndex
    public :: TestAngle2DUtility

    contains

    double precision function Angle2D(K1, K2)
      !Returns the angle in radians between vectors 'K1' and 'K2'
      double precision, dimension(D) :: K1, K2
      double precision :: dotp, det
      !cosang = dot(K1, K2)/norm(K1)/norm(K2)
      dotp=dot_product(K1, K2)
      det=K1(1)*K2(2)-K2(1)*K1(2)
    !   print *, K1, K2
    !   print *, dotp, det
      Angle2D=atan2(det, dotp)
      if(Angle2D<0) Angle2D=Angle2D+2.0*pi
      return
    end function

    double precision function Index2Angle(Index, AngleNum)
      !Map index [1...AngleNum] to the theta range [0.0, 2*pi)
      integer :: Index, AngleNum
      if(Index>AngleNum .or. Index<1) then
        print *, "Angle Index must be between [1, AngleNum]"
        stop
      endif
      Index2Angle=(Index-1)*2.0*pi/AngleNum
      return
    end function

    integer function Angle2Index(Angle, AngleNum)
      !Map theta range  [0.0, 2*pi) to index [1...AngleNum]
      integer :: AngleNum
      double precision :: Angle, temp, dAngle
      if(Angle>2.0*pi .or. Angle<0.0) then
        print *, "Angle must be between [0.0, 2.0*pi)"
        stop
      endif
      dAngle=2.0*pi/AngleNum
      if(Angle>=2.0*pi-dAngle/2.0 .or. Angle<dAngle/2.0) then
        Angle2Index=1
      else
        temp=Angle/dAngle+0.5
        Angle2Index=int(temp)+1
      endif
      return
    end

    subroutine TestAngle2DUtility()
        implicit none
        call TestAngle2D()
        call TestAngleIndex()
    end subroutine

    subroutine TestAngle2D()
        implicit none
        !!!!!!   Test Angle functions  !!!!!!!!!!!!!!!!!!!!!!!
        double precision, dimension(D) :: K1, K2
        K1=0.0; K2=0.0
        K1(1)=1.0; K2(1)=1.0
        if(abs(Angle2D(K1, K2))>1.e-7) then
            print *, "Angle between K1 and K2 are not zero!"
            stop
        endif

        K1=0.0; K2=0.0
        K1(1)=1.0; K2(1)=-1.0
        if(abs(Angle2D(K1, K2)-pi)>1.e-7) then
            print *, "Angle between K1 and K2 are not Pi! Instead, it is ", Angle2D(K1, K2)
            stop
        endif

        K1=0.0; K2=0.0
        K1(1)=1.0; K2(1)=1.0; K2(2)=-eps
        if(abs(Angle2D(K1, K2)-2.0*pi)>1.e-7) then
            print *, "Angle between K1 and K2 are not 2.0*Pi! "
            print *, Angle2D(K1, K2)
            stop
        endif
    end subroutine

    subroutine TestAngleIndex()
        implicit none
        !!!!!!   Test Angle functions  !!!!!!!!!!!!!!!!!!!!!!!
        integer  :: AngleNum=64
        if(abs(Index2Angle(1, AngleNum)-0.0)>eps) then
            print *, "Angle for index 1 should be zero!"
            stop
        endif
        if(abs(Index2Angle(AngleNum, AngleNum)-(2.0*pi*(1.0-1.0/AngleNum)))>1.0e-10) then
            print *, "Angle for index AngleNum should be 2.0*pi-0^+!"
            stop
        endif
        if(Angle2Index(0.d0, AngleNum)/=1) then
            print *, "Angle zero should have index 1!"
            stop
        endif
        if(Angle2Index(2.0*pi*(1.0-0.5/AngleNum)+eps, AngleNum)/=1) then
            print *, "Angle 2*pi-pi/AngleNum should have index 1!"
            stop
        endif
        if(Angle2Index(2.0*pi*(1.0-0.5/AngleNum)-eps, AngleNum)/=AngleNum) then
            print *, "Angle 2*pi-pi/AngleNum-0^+ should have index AngleNum!"
            print *, Angle2Index(2.0*pi*(1.0-0.5/AngleNum)-eps, AngleNum)
            stop
        endif
    end subroutine

end module