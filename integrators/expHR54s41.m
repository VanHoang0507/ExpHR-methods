function [time,expHR54s4_sol,h,CPU,acc,rej] = expHR54s41(F,A,J,t0,t_end,u0,Atol,Rtol,tol,NTS)    
m = length(u0);
zero=zeros(m,1);
time=[];
dt=(t_end-t0)/NTS;
pow = 1/5;
c2=0.5;
c3=3/10;
c4=(c3*5-3)/(10*c3-5);
time(1)=t0;
i=1;
h(1)=dt;
acc=0;
rej=0;
t=t0;
tic

   while t < t_end
        % Ensure the step size does not overshoot the end point
        if t + dt > t_end
            dt = t_end - t;
            time=[time,time(end)+dt];
        end       
        Ju0 = J(u0);
        Fu0 = F(u0);
        Jn=A+Ju0;
        Un2= u0 + phipm_simul_iom(c2*dt,Jn,[zero,Fu0],tol,1,2);
        Hn2=(J(Un2)-Ju0)*F(Un2);
        % U1 = phipm_simul_iom([c3,c4]*dt,Jn,[zero,F(u0),zero,(1/c2)*Hn2/dt],tol,1,2);
        % Un3=u0+U1(:,1);
        % Hn3=(J(Un3)-J(u0))*F(Un3);
        % Un4=u0+U1(:,2);
        % Hn4=(J(Un4)-J(u0))*F(Un4);
        c = [c3,c4];
        Un = zeros(m,2);
        Hn = zeros(m,2);
        maxNumCompThreads(1);
        parfor j=1:2
            Un(:,j) = u0+phipm_simul_iom(c(j)*dt,Jn,[zero,Fu0,zero,1/c2*Hn2/dt],tol,1,2);
            F1 = F(Un(:,j));
            Hn(:,j)=(J(Un(:,j))*F1-Ju0*F1);
        end
        Hn3 = Hn(:,1);
        Hn4 = Hn(:,2);
        u_new=u0+dt^2*phipm_simul_iom(1,dt*Jn,[zero,Fu0/dt,zero,c4*Hn3/(c3*c4-c3^2)+c3*Hn4/(c4*c3-c4^2),-2*Hn3/(c4*c3-c3^2)-2*Hn4/(c4*c3-c4^2)],tol,1,2);
        sc=Atol+max(abs(u0),abs(u_new))*Rtol;
        Er=dt^2*phipm_simul_iom(1,dt*Jn,[zero,zero,zero,zero,(16)*Hn2+(-400 /27)*Hn3+(-128/27)*Hn4],tol,1,2);
        % Estimate the error (for simplicity, use the difference between old and new value)
        error = norm(Er./sc,2)/sqrt(m);        
        % Check if the error is within the tolerance
        if error <= 1
            % Accept the step: update time and solution
            time(i+1) = time(i) + dt;
            h(i+1) = dt;
            t=t+dt;
            i=i+1;
            u0 = u_new;
            acc=acc+1;
        else
            rej=rej+1;
        end
        
        % Adjust the step size based on the error
        if error < 1e-16
            % If error is zero, increase step size aggressively
            s = 2;
        else
            % Compute the scaling factor to adjust step size
            s = 0.8 * (1 / error) ^ pow;
        end
        
        % Limit step size change factor to avoid drastic changes (between 0.1x and 1x)
        dt = dt * min(1.5, max(0.5, s));
    end
CPU=toc;
expHR54s4_sol=u0;
end

     