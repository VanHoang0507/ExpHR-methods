function [time,expRB42s2_sol,h,CPU,acc,rej]=expRB42s2(F,A,g,Jg,t0,t_end,u0,Atol,Rtol,tol,NTS)
m = length(u0);
zero=zeros(m,1);
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
        % Perform a step with exponential Two derivative Rosenbrock
        Jn =@(u) A*u+Jg(u0,u);
        gn = @(u) g(u)-Jg(u0,u);
        U=phipm_simul_iom([3/4,1]*dt,Jn,[zero,F(u0)],tol,1,2);
        Un2=u0+U(:,1);
        Dn2=F(Un2)-Jn(Un2)-gn(u0);
        u_new=u0+phipm_simul_iom(dt,Jn,[zero,F(u0),zero,(32/9*Dn2)/dt^2],tol,1,2); 
        sc=Atol+max(abs(u0),abs(u_new))*Rtol;
        % Estimate the error (for simplicity, use the difference between old and new value)
        error = norm((u_new-u0-U(:,2))./sc,2)/sqrt(m);

        % Check if the error is within the tolerance
        if error <= 1
            % Accept the step: update time and solution
            time(i+1) = time(i) + dt;
            t=t+dt;
            h(i+1) = dt;
            u0 = u_new;
            i=i+1;
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

        % Limit step size change factor to avoid drastic changes (between 0.1x and 2x)
        dt = dt * min(1.5, max(0.5, s));
    end
CPU=toc;
expRB42s2_sol=u0;
end


