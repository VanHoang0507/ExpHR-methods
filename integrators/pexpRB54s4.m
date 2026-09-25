function [time,pexpRB54s4_sol,h,CPU,acc,rej]=pexpRB54s4(F,A,g,J,t0,t_end,u0,Atol,Rtol,tol,NTS)    
m = length(u0);
zero=zeros(m,1);
time=[];
dt=(t_end-t0)/NTS;
pow = 1/5;
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
        % Perform a step with exprb54
      Fu0 = F(u0);
      Jn = @(u) A*u+J(u0,u);
      gn = @(u) g(u)-J(u0,u);
      Un2 = u0+phipm_simul_iom((1/4)*dt,Jn,[zero,Fu0],tol,1,2);
      Dn2 = F(Un2)-Jn(Un2)-gn(u0);
      U1 = phipm_simul_iom([1/2,9/10]*dt,Jn,[zero,F(u0),zero,32*Dn2/dt^2],tol,1,2);
      Un3=u0+U1(:,1);
      Dn3=F(Un3)-Jn(Un3)-gn(u0);
      Un4=u0+U1(:,2);
      Dn4=F(Un4)-Jn(Un4)-gn(u0);
      % c = [1/2,9/10];
      % Un = zeros(m,2);
      % Dn = zeros(m,2);
      % 
      % parfor j=1:2
      %     Un(:,j) = u0+phipm_simul_iom(c(j)*dt,Jn,[zero,Fu0,zero,32*Dn2/dt^2],tol,1,2);
      %     Dn(:,j) = F(Un(:,j))-Jn(Un(:,j))-gn(u0);
      % end
      % Dn3 = Dn(:,1);
      % Dn4 = Dn(:,2);
      
      u_new=u0+phipm_simul_iom(dt,Jn,[zero,Fu0,zero,(18*Dn3-(250/81)*Dn4)/dt^2,(-60*Dn3+(500/27)*Dn4)/dt^3],tol,1,2);
      sc = Atol+max(abs(u0),abs(u_new))*Rtol;
      Er = phipm_simul_iom(dt,Jn,[zero,zero,zero,(-64*Dn2+26*Dn3-(250/81)*Dn4)/dt^2,(60*Dn2-(195/8)*Dn3+(625/216)*Dn4)/dt^3],tol,1,2);
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
        
        % Limit step size change factor to avoid drastic changes (between 0.1x and 2x)
        dt = dt * min(1.5, max(0.5, s));
    end
CPU=toc;
pexpRB54s4_sol=u0;
end

     