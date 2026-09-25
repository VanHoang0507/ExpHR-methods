function [time,pexpRB54s5_sol,h,CPU,acc,rej]=pexpRB54s5(F,A,g,J,t0,t_end,u0,Atol,Rtol,tol,NTS)    
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
        % Perform a step with exprb55
      Jn =@(u) A*u+J(u0,u);
      gn=@(u) g(u)-J(u0,u);
      Fu0 = F(u0);
      Un2=u0+phipm_simul_iom(1/2*dt,Jn,[zero,Fu0],tol,1,2);
      Dn2=F(Un2)-Jn(Un2)-gn(u0);
      U1 = phipm_simul_iom([1/3,1/2,1]*dt,Jn,[zero,Fu0,zero,8*Dn2/dt^2],tol,1,2);
      Un3=u0+U1(:,2);
      Dn3=F(Un3)-Jn(Un3)-gn(u0);
      Un4=u0+U1(:,1);
      Dn4=F(Un4)-Jn(Un4)-gn(u0);
      Un5=u0+U1(:,3);
      Dn5=F(Un5)-Jn(Un5)-gn(u0);
      % c = [1/3,1/2,1];
      % Un = zeros(m,3);
      % Dn = zeros(m,3);
      % maxNumCompThreads(1);
      % parfor j=1:3
      %     Un(:,j) = u0+phipm_simul_iom(c(j)*dt,Jn,[zero,Fu0,zero,8*Dn2/dt^2],tol,1,2);
      %     Dn(:,j)=F(Un(:,j))-Jn(Un(:,j))-gn(u0);
      % end
      % Dn3 = Dn(:,2);
      % Dn4 = Dn(:,1);
      % Dn5 = Dn(:,3);
      u_new=u0+phipm_simul_iom(dt,Jn,[zero,Fu0,zero,(-32*Dn3+81*Dn4+Dn5)/dt^2,(384*Dn3-729*Dn4-15*Dn5)/dt^3,(-1152*Dn3+1944*Dn4+72*Dn5)/dt^4],tol,1,2);
      sc=Atol+max(abs(u0),abs(u_new))*Rtol;
      Er=phipm_simul_iom(dt,Jn,[zero,zero,zero,zero,zero,(-1152*Dn3+1944*Dn4+72*Dn5)/dt^4],tol,1,2);
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
pexpRB54s5_sol=u0;
end

     