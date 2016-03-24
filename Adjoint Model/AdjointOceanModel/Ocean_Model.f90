    module Ocean_Model
    use Const
    implicit none
    
    ! 正向积分模型中的水位变化
    real(kind=PREC), allocatable, &
    dimension(:, :, :)                  :: H, U, V
    
    ! 用于调和分析累加求和的临时变量
    ! 以及算得的整个模拟区域的调和常数
    real(kind=PREC), allocatable, &
    dimension(:, :)                     :: tmp_Tide_A, tmp_Tide_B, simZeta, simSigma
    
    contains
    
    subroutine Ocean_Model_allocate
        implicit none
        
        allocate(H(1:numSteps,numX,numY))
        allocate(U(-2:numSteps,numX,numY))
        allocate(V(-2:numSteps,numX,numY))
        allocate(tmp_Tide_A(numX,numY))
        allocate(tmp_Tide_B(numX,numY))
        allocate(simZeta(numX,numY))
        allocate(simSigma(numX,numY))
        
    end subroutine
    
    subroutine Ocean_Model_Run
        implicit none
        real(kind=PREC) :: maxDiff, temp(56, 43)
        integer         :: iter, t
        
        ! 初始化
        H = 0.
        U = 0.
        V = 0.
        
        do iter = 1, maxIteration
            
            ! 首先把上一次模拟的末状态复制到当前的开始状态
            U(1, :, :) = U(numSteps, :, :)
            V(1, :, :) = V(numSteps, :, :)
            H(1, :, :) = H(numSteps, :, :)
            
            do t = 1, numSteps - 2, 2
                ! 交错步迭代
                call calc_H(t)
                call calc_V(t, t)
                call calc_U(t, t + 1)
                call calc_H(t + 1)
                call calc_U(t + 1, t + 1)
                call calc_V(t + 1, t + 2)
            end do
            
            ! 计算最大误差
            maxDiff = maxval(dabs(H(numSteps, :, :) - H(1, :, :)))
            write(*, "(4X, 'Error = ', E12.5E2, ' on forward iteration:', I3)") maxDiff, iter
            ! 检查是否满足退出条件
            if ( maxDiff < IntegrationEps ) exit
            
        end do
        
    end subroutine
    
    subroutine calc_H( curStep )
        implicit none
        integer         :: curStep
        real(kind=PREC) :: depthL, depthR, depthU, depthD
        
        do x=2,numX-1  
            do y=2,numY-1
                if (maskBoundary(x,y) == 1) then
                        depthL=0.5*(WaterDepth(x+1,y)+WaterDepth(x,y))
                        depthR=0.5*(WaterDepth(x,y)+WaterDepth(x-1,y))
                        depthU=0.5*(WaterDepth(x,y+1)+WaterDepth(x,y))
                        depthD=0.5*(WaterDepth(x,y)+WaterDepth(x,y-1))
                        H(curStep+1,x,y)=H(curStep,x,y) &
                                        -ddx(y)*(depthL*U(curStep,x,y)-depthR*U(curStep,x-1,y)) &
                                        -ddy*(depthU*V(curStep,x,y)-depthD*V(curStep,x,y-1))
                    else if (maskBoundary(x,y) > 1) then
                        H(curStep+1,x,y)=Tide_A(x,y)*dcos(TideOmega*curStep*dt) &
                                        +Tide_B(x,y)*dsin(TideOmega*curStep*dt)
                    else 
                        H(curStep+1,x,y)=0.
                    end if
            end do
        end do
    end subroutine
    
    subroutine calc_U( curStep, relyStep )
        implicit none
        integer         :: curStep, relyStep
        real(kind=PREC) :: avq, ap0, ap1, r0, ee, ff, aa, bb, cc, dd
        
        do y=2,numY-1
            do x=2,numX-1
                if (maskU(x,y) == 1) then
                    avq=maskV(x,y)+maskV(x+1,y)+maskV(x,y-1)+maskV(x+1,y-1)
                    ap0=(V(curStep,x,y)+V(curStep,x+1,y)+V(curStep,x,y-1)+V(curStep,x+1,y-1))/max(1.0,avq)
                    ap1=(V(relyStep,x,y)+V(relyStep,x+1,y)+V(relyStep,x,y-1)+V(relyStep,x+1,y-1))/max(1.0,avq)
                    r0=dsqrt(U(curStep,x,y)**2+ap0**2)  
                    ee=0.5*(WaterDepth(x,y)+WaterDepth(x+1,y))
                    ff=0.5*(Friction(x,y)+Friction(x+1,y))

                    aa=1.+alpha*dt*ff*r0/ee
                    bb=1.-(1.-alpha)*dt*ff*r0/ee
                    cc=(U(curStep,x+1,y)-2.*U(curStep,x,y)+U(curStep,x-1,y))/dx(y)**2 &
                    +(U(curStep,x,y+1)-2.*U(curStep,x,y)+U(curStep,x,y-1))/dy**2
                    dd=U(curStep,x,y)*(U(curStep,x+1,y)-U(curStep,x-1,y))/dx(y)/2. &
                    +ap1*(U(curStep,x,y+1)-U(curStep,x,y-1))/dy/2.

                    U(curStep+1,x,y)=(bb/aa)*U(curStep,x,y) &
                    -(ddx(y)/aa)*g*(H(curStep+1,x+1,y)-H(curStep+1,x,y)) &
                    +(dt/aa)*Viscosity(x,y)*cc &
                    +(dt/aa)*f(y)*ap1 &
                    -(dt/aa)*dd
                else
                    U(curStep+1,x,y)=0.
                end if
            end do
        end do
    end subroutine
    
    subroutine calc_V( curStep, relyStep )
        implicit none
        integer         :: curStep, relyStep
        real(kind=PREC) :: avp, aq0, aq1, s0, ee, ff, aa, bb, cc, dd
        
        do x=2,numX-1
            do y=2,numY-1
                if (maskV(x,y) == 1) then
                    avp=maskU(x,y)+maskU(x-1,y)+maskU(x,y+1)+maskU(x-1,y+1)
                    aq0=(U(curStep,x,y)+U(curStep,x-1,y)+U(curStep,x,y+1)+U(curStep,x-1,y+1))/max(1.0,avp)
                    aq1=(U(relyStep,x,y)+U(relyStep,x-1,y)+U(relyStep,x,y+1)+U(relyStep,x-1,y+1))/max(1.0,avp)
                    s0=dsqrt(aq0**2+V(curStep,x,y)**2)  
                    ee=0.5*(WaterDepth(x,y)+WaterDepth(x,y+1))
                    ff=0.5*(Friction(x,y)+Friction(x,y+1)) 

                    aa=1.+alpha*dt*ff*s0/ee
                    bb=1.-(1.-alpha)*dt*ff*s0/ee
                    cc=(V(curStep,x+1,y)-2.*V(curStep,x,y)+V(curStep,x-1,y))/dx(y)**2 &
                    +(V(curStep,x,y+1)-2.*V(curStep,x,y)+V(curStep,x,y-1))/dy**2
                    dd=aq1*(V(curStep,x+1,y)-V(curStep,x-1,y))/dx(y)/2. &
                    +V(curStep,x,y)*(V(curStep,x,y+1)-V(curStep,x,y-1))/dy/2.

                    V(curStep+1,x,y)=(bb/aa)*V(curStep,x,y) &
                    -(ddy/aa)*g*(H(curStep+1,x,y+1)-H(curStep+1,x,y)) &
                    +(dt/aa)*Viscosity(x,y)*cc &
                    -(dt/aa)*f(y)*aq1 &
                    -(dt/aa)*dd
                else
                    V(curStep+1,x,y)=0.
                end if
            end do
        end do
    end subroutine
    
    ! 调和分析：根据模拟结果，计算整个模拟区域的调和常数
    ! 然后与观测点处的调和常数对比
    subroutine Harmonic
        implicit none
        integer :: t
        
        tmp_Tide_A = 0.
        tmp_Tide_B = 0.
        do t = 2, numSteps
            tmp_Tide_A = tmp_Tide_A + H(t, :, :) * dcos(TideOmega*(t-1)*dt)
            tmp_Tide_B = tmp_Tide_B + H(t, :, :) * dsin(TideOmega*(t-1)*dt)
        end do
        tmp_Tide_A = 2. * tmp_Tide_A / real((numSteps - 1))
        tmp_Tide_B = 2. * tmp_Tide_B / real((numSteps - 1))
        
        do x = 1, numX
            do y = 1, numY
                simZeta(x, y) = dsqrt(tmp_Tide_A(x, y)**2 + tmp_Tide_B(x, y)**2)
                if ( simZeta(x, y) < eps ) then
                    simSigma(x, y) = 0.
                else
                    if ( tmp_Tide_B(x, y) > 0 ) then
                        simSigma(x, y) = dacos(tmp_Tide_A(x, y) / simZeta(x, y))
                    else
                        simSigma(x, y) = 2. * PI - dacos(tmp_Tide_A(x, y) / simZeta(x, y))
                    end if
                    simSigma(x, y) = simSigma(x, y) * 180. / PI
                end if
            end do
        end do
        
    end subroutine

    end module