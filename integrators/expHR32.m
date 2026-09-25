function [t, expHR32_sol,h,CPU,acc,rej] = expHR32(F,A,J,t0,t_end,u0,Atol,Rtol,tol,c,NTS)
m = length(u0);
zero=zeros(m,1);
time=[];
dt=(t_end-t0)/NTS;
pow = 1/3;
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
        % Perform a step with exprb43
        Jn=(A+J(u0));         
        Un=phipm_simul_iom([c,1]*dt,Jn,[zero,F(u0)],tol,1,2);
        Hn2=(J(u0+Un(:,1))-J(u0))*F(Un(:,1));
        u_new=u0+Un(:,2)+dt*phipm_simul_iom(1,dt*Jn,[zero,zero,zero,dt*(1/c)*Hn2],tol,1,2);
        sc=Atol+max(abs(u0),abs(u_new))*Rtol;
        Er=norm((dt^2*phipm_simul_iom(1,dt*Jn,[zero,zero,zero,(1/c)*Hn2],tol,1,2))./sc,2)/sqrt(m);
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
            s = 0.9 * (1 / error) ^ pow;
        end
        
        % Limit step size change factor to avoid drastic changes (between 0.1x and 2x)
        dt = dt * min(2, max(0.1, s));
    end

CPU=toc;
expHR32_sol=u0;
end