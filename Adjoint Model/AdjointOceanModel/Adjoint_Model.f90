    module Adjoint_Model
    use Const
    use Ocean_Model
    implicit none
    
    real(kind=PREC), allocatable, &
    dimension(:)                    :: costFunction
    
    real(kind=PREC), allocatable, &
    dimension(:, :, :)              :: diffH, H_adj, V_adj, U_adj, ap, aq, R, S
    
    real(kind=PREC), allocatable, &
    dimension(:, :)                 :: diffZeta, diffSigma, pJ_pa, pJ_pb, pJ_pf
    
    real(kind=PREC)                 :: avgSigma, avgZeta, rmsSigma, rmsZeta
    
    contains
    
    subroutine Adjoint_Model_allocate
        implicit none
        
        allocate(costFunction(maxIteration))
        allocate(diffH(numSteps, numX, numY))
        allocate(H_adj(numSteps, numX, numY))
        allocate(V_adj(numSteps, numX, numY))
        allocate(U_adj(numSteps, numX, numY))
        allocate(ap(-3:0,numX,numY))
        allocate(aq(-3:0,numX,numY))
        allocate(R(-3:0,numX,numY))
        allocate(S(-3:0,numX,numY))
        allocate(diffZeta(numX,numY))
        allocate(diffSigma(numX,numY))
        allocate(pJ_pa(numX,numY))
        allocate(pJ_pb(numX,numY))
        allocate(pJ_pf(numX,numY))
        
        costFunction = 0.
        
    end subroutine
    
    subroutine Adjoint_Model_Run
        implicit none
        integer         :: iter, t
        
        do iter = 1, maxIteration
            write(*, '("Current optimize iteration:", I3)') iter
            
            ! 应用先前的优化结果
            ! 注意只处理开边界上面的调和常数
            where ( maskBoundary > 1 )
                Tide_A = newTide_A
                Tide_B = newTide_B
            end where
            
            Friction = newFriction
            
            ! 运行正向积分模型
            call Ocean_Model_Run
            
            ! 正向积分完成后，开始计算模拟结果与观测点的差异
            do t = 1, numSteps
                diffH(t, :, :) = maskObserve * ( H(t, :, :) - H0(t, :, :) )
            end do
            !write(*, *) maxval(dabs(diffH))
            !call exit(1)
            ! 然后计算代价函数
            costFunction(iter) = sum(diffH ** 2) * 0.5
            write(*, "('Value of cost function:', F15.6)") costFunction(iter)
            write(*, "('Cost function decrease:', F15.6)") costFunction(iter) / costFunction(1)
            
            ! 如果代价函数出现波动，就自动减小步长
            if ( iter > 1 &
                .and. costFunction(iter) > costFunction(iter - 1) &
                .and. alphaTideAB > alphaTideABmin ) then
                alphaTideAB = alphaTideAB - alphaTideABdec;
                write(*, "('Alpha of a, b is decreased, current value is: ', E12.4)") alphaTideAB
            end if
            
            ! 进行调和分析，获得整个模拟区域的调和常数
            call Harmonic
            
            ! 计算模拟值和观测点值的差异
            call Harmonic_Diff
            
            ! 初始化伴随变量
            H_adj = 0.
            U_adj = 0.
            V_adj = 0.
            
            ! 运行伴随模式
            do t = numSteps, 3, -2
                call Adjoint_init(t)
                call Adjoint_H(t, -1)
                call Adjoint_V(t, t, -1)
                call Adjoint_U(t, t - 1, -1)
                call Adjoint_H(t - 1, -2)
                call Adjoint_U(t - 1, t - 1, -2)
                call Adjoint_V(t - 1, t - 2, -2)
            end do
            
            ! 计算相应参数的梯度并优化
            if ( optOpenBoundary ) call Grad_AB
            if ( optFrictionData ) call Grad_Friction
            
            ! TODO：修改一下输出格式
            ! 输出差异的绝对值平均和方均根
            ! write(*, *) avgSigma, avgZeta
            write(*, "('Root-mean-square of Zeta, Sigma:', 2F15.6)") rmsZeta, rmsSigma
            write(*, "('Average absolute of Zeta, Sigma:', 2F15.6)") avgZeta, avgSigma
            write(*, '("**************************************************************")') 
            
            if ( ExitType == 0 ) then
                if ( rmsSigma < SigmaEps .and. rmsZeta < ZetaEps ) exit
            else if ( ExitType == 1 ) then
                if ( avgSigma < SigmaEps .and. avgZeta < ZetaEps ) exit
            else
                write(*, *) "Error: unknown exit type!"
                call exit(3)
            end if
            
        end do
        ! 输出代价函数的变化
        call dump1D(costFunction, 'costFunction.dat', 1, iter, maxIteration)
        
    end subroutine
    
    subroutine Harmonic_Diff
        implicit none
        
        diffZeta = maskObserve * ( simZeta - obsZeta )
        diffSigma = maskObserve * ( simSigma - obsSigma )
        
        ! 将角度区间修正到-Pi~Pi
        where ( diffSigma > 180. )
            diffSigma = diffSigma - 360.
        else where ( diffSigma < -180. )
            diffSigma = diffSigma + 360.
        end where
        
        ! 对于无潮点，差异为零
        where ( simZeta < eps .and. obsZeta < eps )
            diffSigma = 0.
        end where
        
        ! 求差异平均值和方均根
        avgSigma = sum(dabs(diffSigma)) / real(sumObserve)
        rmsSigma = dsqrt(sum(diffSigma**2) / real(sumObserve))
        avgZeta = sum(dabs(diffZeta)) / real(sumObserve)
        rmsZeta = dsqrt(sum(diffZeta**2) / real(sumObserve))
        
    end subroutine
    
    subroutine Adjoint_init( curStep )
        implicit none
        integer         :: i, curStep
        real(kind=PREC) :: avp, avq
        ap = 0.
        aq = 0.
        R = 0.
        S = 0.
        do x=2,numX-1
            do y=2,numY-1
                avp=maskU(x,y)+maskU(x-1,y)+maskU(x,y+1)+maskU(x-1,y+1)
                do i=1,4
                    ap(1-i,x,y)=(U(curStep+1-i,x,y)+U(curStep+1-i,x-1,y) &
                                +U(curStep+1-i,x,y+1)+U(curStep+1-i,x-1,y+1))/max(1.0,avp)
                    S(1-i,x,y)=dsqrt(ap(1-i,x,y)**2+V(curStep+1-i,x,y)**2)
                end do
                avq=maskV(x,y)+maskV(x+1,y)+maskV(x,y-1)+maskV(x+1,y-1)
                do i=1,4
                    aq(1-i,x,y)=(V(curStep+1-i,x,y)+V(curStep+1-i,x+1,y) &
                                +V(curStep+1-i,x,y-1)+V(curStep+1-i,x+1,y-1))/max(1.0,avq)
                    R(1-i,x,y)=dsqrt(aq(1-i,x,y)**2+U(curStep+1-i,x,y)**2)
                end do
            end do
        end do
    end subroutine
    
    subroutine Adjoint_H( curStep, extStep )
        implicit none
        integer         :: curStep, extStep
        real(kind=PREC) :: ff1, ff2, ee1, ee2
        
        do x=2,numX-1
            do y=2,numY-1
                if (maskBoundary(x,y) == 1) then
                    ff1=0.5*(Friction(x,y)+Friction(x+1,y))
                    ff2=0.5*(Friction(x,y)+Friction(x,y+1))
                    ee1=0.5*(WaterDepth(x,y)+WaterDepth(x+1,y))
                    ee2=0.5*(WaterDepth(x,y)+WaterDepth(x,y+1))
                    H_adj(curStep-1,x,y)=H_adj(curStep,x,y) &
                                        -ddx(y)*g*(U_adj(curStep,x-1,y)-U_adj(curStep,x,y)) &
                                        -ddy*g*(V_adj(curStep,x,y-1)-V_adj(curStep,x,y)) &
                                        -ddx(y)*U(curStep-1,x,y)*(H_adj(curStep,x-1,y)-H_adj(curStep,x,y)) &
                                        -ddy*V(curStep-1,x,y)*(H_adj(curStep,x,y-1)-H_adj(curStep,x,y)) &
                                        +ff1*U_adj(curStep,x,y)*U(curStep-1,x,y)*R(extStep,x,y)/ee1**2 &
                                        +ff2*V_adj(curStep,x,y)*V(curStep-1,x,y)*S(extStep,x,y)/ee2**2 &
                                        -maskObserve(x,y)*diffH(curStep,x,y)
                else
                    H_adj(curStep-1,x,y)=0.
                end if
            end do
        end do
        
    end subroutine
    
    subroutine Adjoint_U( curStep, relyStep, extStep )
        implicit none
        integer         :: curStep, relyStep, extStep
        real(kind=PREC) :: avq, aqa_1, eea, ffa, gga, dda, aaa, bba, cca, ppq
        do x=2,numX-1
            do y=2,numY-1
                if (maskU(x,y) == 1) then
                    avq=maskV(x,y)+maskV(x+1,y)+maskV(x,y-1)+maskV(x+1,y-1)
                    aqa_1=(V_adj(relyStep,x,y)+V_adj(relyStep,x+1,y) &
                        +V_adj(relyStep,x,y-1)+V_adj(relyStep,x+1,y-1))/max(1.0,avq)
                    eea=0.5*(WaterDepth(x,y)+WaterDepth(x+1,y))
                    ffa=0.5*(Friction(x,y)+Friction(x+1,y))
                    gga=(U_adj(curStep,x+1,y)-2.*U_adj(curStep,x,y)+U_adj(curStep,x-1,y))/dx(y)**2 &
                        +(U_adj(curStep,x,y+1)-2.*U_adj(curStep,x,y)+U_adj(curStep,x,y-1))/dy**2
                    dda=aqa_1*(V(curStep-1,x+1,y)-V(curStep-1,x-1,y))/dx(y)/2. &
                        +U_adj(curStep,x,y)*(U(curStep-1,x+1,y)-U(curStep-1,x-1,y))/dx(y)/2. &
                        -(U_adj(curStep,x,y+1)*aq(extStep,x,y+1)-U_adj(curStep,x,y-1)*aq(extStep,x,y-1))/dy/2. &
                        -(U_adj(curStep,x+1,y)*U(curStep-1,x+1,y)-U_adj(curStep,x-1,y)*U(curStep-1,x-1,y))/dx(y)/2.

                    if (R(extStep,x,y) < eps) then
                        aaa=1.
                        bba=1.
                        cca=f(y)
                    else
                        ppq=ffa*(aq(extStep,x,y)**2+2.*U(curStep-1,x,y)**2)/R(extStep,x,y)/eea
                        aaa=1+(1-alpha)*dt*ppq
                        bba=1-alpha*dt*ppq
                        cca=f(y)+ffa*aq(extStep,x,y)*U(curStep-1,x,y)/R(extStep,x,y)/eea
                    end if

                    U_adj(curStep-1,x,y)=(bba/aaa)*U_adj(curStep,x,y) &
                                        +(ddx(y)/aaa)*eea*(H_adj(curStep-1,x+1,y)-H_adj(curStep-1,x,y)) &
                                        -(dt/aaa)*cca*aqa_1 &
                                        +(dt/aaa)*Viscosity(x,y)*gga &
                                        -(dt/aaa)*dda  
                else
                    U_adj(curStep-1,x,y)=0.
                end if
            end do
        end do
    end subroutine
    
    subroutine Adjoint_V( curStep, relyStep, extStep )
        implicit none
        integer         :: curStep, relyStep, extStep
        real(kind=PREC) :: avp, apa_1, eea, ffa, gga, dda, aaa, bba, cca, pqq
        
        do x=2,numX-1
            do y=2,numY-1
                if (maskV(x,y) == 1) then
                    avp=maskU(x,y)+maskU(x-1,y)+maskU(x,y+1)+maskU(x-1,y+1)
                    apa_1=(U_adj(relyStep,x,y)+U_adj(relyStep,x-1,y) &
                        +U_adj(relyStep,x,y+1)+U_adj(relyStep,x-1,y+1))/max(1.0,avp)
                    eea=0.5*(WaterDepth(x,y)+WaterDepth(x,y+1))
                    ffa=0.5*(Friction(x,y)+Friction(x,y+1))
                    gga=(V_adj(curStep,x+1,y)-2.*V_adj(curStep,x,y)+V_adj(curStep,x-1,y))/dx(y)**2 &
                        +(V_adj(curStep,x,y+1)-2.*V_adj(curStep,x,y)+V_adj(curStep,x,y-1))/dy**2
                    dda=apa_1*(U(curStep-1,x,y+1)-U(curStep-1,x,y-1))/dy/2. &
                        +V_adj(curStep,x,y)*(V(curStep-1,x,y+1)-V(curStep-1,x,y-1))/dy/2. &
                        -(V_adj(curStep,x+1,y)*ap(extStep,x+1,y)-V_adj(curStep,x-1,y)*ap(extStep,x-1,y))/dx(y)/2. &
                        -(V_adj(curStep,x,y+1)*V(curStep-1,x,y+1)-V_adj(curStep,x,y-1)*V(curStep-1,x,y-1))/dy/2.

                    if (S(extStep,x,y) < eps) then
                        aaa=1.
                        bba=1.
                        cca=f(y)
                    else
                        pqq=ffa*(ap(extStep,x,y)**2+2.*V(curStep-1,x,y)**2)/S(extStep,x,y)/eea
                        aaa=1.+(1.-alpha)*dt*pqq
                        bba=1.-alpha*dt*pqq
                        cca=f(y)-ffa*ap(extStep,x,y)*V(curStep-1,x,y)/S(extStep,x,y)/eea
                    end if

                    V_adj(curStep-1,x,y)=(bba/aaa)*V_adj(curStep,x,y) &
                                        +(ddy/aaa)*eea*(H_adj(curStep-1,x,y+1)-H_adj(curStep-1,x,y)) &
                                        +(dt/aaa)*cca*apa_1 &
                                        +(dt/aaa)*Viscosity(x,y)*gga &
                                        -(dt/aaa)*dda  	  
                else
                    V_adj(curStep-1,x,y)=0.
                end if
            end do
        end do
        
    end subroutine
    
    subroutine Grad_Friction
        implicit none
        real(kind=PREC)                         :: avp, avq, r_1, s_1, ap_1, aq_1, sum_F
        real(kind=PREC), dimension(numX, numY)  :: PJKP, PJKQ
        
        PJKP = 0.
        PJKQ = 0.
        
        do x=2,numX-1
            do y=2,numY-1
                if (maskU(x,y) == 1 .and. maskV(x,y) == 1) then
                    do t=2,numSteps-1
                        avq=maskV(x,y)+maskV(x+1,y)+maskV(x,y-1)+maskV(x+1,y-1)
                        aq_1=(V(t-1,x,y)+V(t-1,x+1,y) &
                            +V(t-1,x,y-1)+V(t-1,x+1,y-1))/max(1.0,avq)
                        r_1=dsqrt(aq_1**2+U(t-1,x,y)**2)
                        PJKP(x,y)=PJKP(x,y)+U_adj(t,x,y)*r_1 &
                                    *((1-alpha)*U(t,x,y)+alpha*U(t+1,x,y)) &
                                    /(0.5*(WaterDepth(x+1,y)+WaterDepth(x,y)))
                    end do
                end if
            end do
        end do

        do x=2,numX-1
            do y=2,numY-1
                if (maskU(x,y) == 1 .and. maskV(x,y) == 1) then
                    do t=2,numSteps-1
                        avp=maskU(x,y)+maskU(x-1,y)+maskU(x,y+1)+maskU(x-1,y+1)
                        ap_1=(U(t-1,x,y)+U(t-1,x-1,y)&
                            +U(t-1,x,y+1)+U(t-1,x-1,y+1))/max(1.0,avp)
                        s_1=dsqrt(ap_1**2+V(t-1,x,y)**2)
                        PJKQ(x,y)=PJKQ(x,y)+V_adj(t,x,y)*s_1 &
                                    *((1-alpha)*V(t,x,y)+alpha*V(t+1,x,y)) &
                                    /(0.5*(WaterDepth(x,y+1)+WaterDepth(x,y)))
                    end do  
                end if	  
            end do
        end do
        
        pJ_pf = PJKP + PJKQ   
        sum_F = sum(pJ_pf ** 2)
        pJ_pf = pJ_pf / dsqrt(sum_F)
        newFriction = Friction - alphaFriction * pJ_pf
        ! 保证底摩擦是正数
        newFriction = maskH * dabs(newFriction)
        write(*, "('Sum of gradient friction:', F15.6)") sum_F
    end subroutine
    
    subroutine Grad_AB
        implicit none
        real(kind=PREC) :: sum_AB
        pJ_pa = 0.
        pJ_pb = 0.
        do t=1,numSteps-1
            do y=2,numY-1
                do x=2,numX-1
                    if (maskBoundary(x,y) == 2) then
                        pJ_pa(x,y)=pJ_pa(x,y)-g*U_adj(t,x,y)*dcos(TideOmega*t*dt)/dx(y)
                        pJ_pb(x,y)=pJ_pb(x,y)-g*U_adj(t,x,y)*dsin(TideOmega*t*dt)/dx(y)
                    else if (maskBoundary(x,y) == 3) then
                        pJ_pa(x,y)=pJ_pa(x,y)+g*V_adj(t,x,y-1)*dcos(TideOmega*t*dt)/dy
                        pJ_pb(x,y)=pJ_pb(x,y)+g*V_adj(t,x,y-1)*dsin(TideOmega*t*dt)/dy
                    else if (maskBoundary(x,y) == 4) then
                        pJ_pa(x,y)=pJ_pa(x,y)+g*U_adj(t,x-1,y)*dcos(TideOmega*t*dt)/dx(y)
                        pJ_pb(x,y)=pJ_pb(x,y)+g*U_adj(t,x-1,y)*dsin(TideOmega*t*dt)/dx(y)
                    else if (maskBoundary(x,y) == 5) then
                        pJ_pa(x,y)=pJ_pa(x,y)-g*V_adj(t,x,y)*dcos(TideOmega*t*dt)/dy
                        pJ_pb(x,y)=pJ_pb(x,y)-g*V_adj(t,x,y)*dsin(TideOmega*t*dt)/dy
                    else if (maskBoundary(x,y) == 6) then
                        pJ_pa(x,y)=pJ_pa(x,y)-g*V_adj(t,x,y)*dcos(TideOmega*t*dt)/dy &
                                    +g*U_adj(t,x-1,y)*dcos(TideOmega*t*dt)/dx(y)
                        pJ_pb(x,y)=pJ_pb(x,y)-g*V_adj(t,x,y)*dsin(TideOmega*t*dt)/dy &
                                    +g*U_adj(t,x-1,y)*dsin(TideOmega*t*dt)/dx(y)
                    end if
                end do
            end do
        end do
        sum_AB = sum(pJ_pa ** 2 + pJ_pb ** 2)
        pJ_pa = pJ_pa / dsqrt(sum_AB)
        pJ_pb = pJ_pb / dsqrt(sum_AB)
        newTide_A = Tide_A - alphaTideAB * pJ_pa
        newTide_B = Tide_B - alphaTideAB * pJ_pb
        write(*, "('Sum of gradient A and B:', F15.6)") sum_AB
        
    end subroutine
    
    end module