function [time,expHR43s3_sol,h,CPU,acc,rej] = expHR43s3(F,A,J,t0,t_end,u0,Atol,Rtol,tol,NTS)    
m = length(u0);
zero=zeros(m,1);
time=[];
dt=(t_end-t0)/NTS;
pow = 1/4;
time(1)=t0;
i=1;
h(1)=dt;
acc=0;
rej=0;
t=t0;
c2 = 0.1;
c3 = 0.6;
a32 =  c3/(c2*(c3-c2));
a33 =  c2/(c3*(c2-c3));

a42 = -2/(c2*(c3-c2));
a43 = -2/(c3*(c2-c3));
tic
   while t < t_end
        % Ensure the step size does not overshoot the end point
        if t + dt > t_end
            dt = t_end - t;
            time=[time,time(end)+dt];
        end       
        % Perform a step with exprb43
        Jn = @(u)(A*u+J(u0,u));
        U=phipm_simul_iom([c2,c3]*dt,Jn,[zero,F(u0)],tol,1,2);
        Un2=u0+U(:,1);
        FUn2 = F(Un2);
        Hn2=(J(Un2,FUn2)-J(u0,FUn2));
        Un3=u0+U(:,2);
        FUn3 = F(Un3);
        Hn3=(J(Un3,FUn3)-J(u0,FUn3));
        phi3_coeff = a32*Hn2 + a33*Hn3;
        phi4_coeff = a42*Hn2 + a43*Hn3;

        u_new=u0+phipm_simul_iom(dt,Jn,[zero,F(u0),zero,phi3_coeff/dt,phi4_coeff/dt^2],tol,1,2);
        sc=Atol+max(abs(u0),abs(u_new))*Rtol;
        Er = phipm_simul_iom(dt,Jn,[zero,zero,zero,zero,phi4_coeff/dt^2],tol,1,2);
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
            s = 1.5;
        else
            % Compute the scaling factor to adjust step size
            s = 0.8 * (1 / error) ^ pow;
        end
        
        % Limit step size change factor to avoid drastic changes (between 0.1x and 2x)
        dt = dt * min(1.5, max(0.5, s));
    end
CPU=toc;
expHR43s3_sol=u0;
end

     