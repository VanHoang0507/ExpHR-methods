function [time,expHR54s4a_sol,h,CPU,acc,rej] = expHR54s4a(F,A,J,t0,t_end,u0,Atol,Rtol,tol,NTS)    
m = length(u0);
zero=zeros(m,1);
time=[];
dt=(t_end-t0)/NTS;
pow = 1/5;
c2 = 0.1;
c3 = 0.7;
c4 = (c3*5-3)/(10*c3-5);
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
        Fu0 = F(u0);
        Jn =@(u) A*u+J(u0,u);
       
        
        U = phipm_simul_iom([c2,c4,c3]*dt,Jn,[zero,Fu0],tol,1,2);
        
        Un2 = u0 + U(:,1);
        FUn2 = F(Un2);
        Hn2=(J(Un2,FUn2)-J(u0,FUn2));
        
        Un3 = u0 + U(:,3) + dt^2*(c3^3/c2)*Hn2/6 ;
        FUn3 = F(Un3);
        Hn3 = (J(Un3,FUn3)-J(u0,FUn3));
        Un4 = u0 + U(:,2) + dt^2*(c4^3/c2)*Hn2/6;
        FUn4 = F(Un4);
        Hn4=(J(Un4,FUn4)-J(u0,FUn4));
        % c = [c3,c4];
        % Un = zeros(m,2);
        % Hn = zeros(m,2);
        % maxNumCompThreads(1);
        % parfor j=1:2
        %     Un(:,j) = u0+phipm_simul_iom(c(j)*dt,Jn,[zero,Fu0,zero,1/c2*Hn2/dt],tol,1,2);
        %     F1 = F(Un(:,j));
        %     Hn(:,j)=J(Un(:,j),F1)-J(u0,F1);
        % end
        % Hn3 = Hn(:,1);
        % Hn4 = Hn(:,2);
        u_new=u0+phipm_simul_iom(dt,Jn,[zero,Fu0,zero,(c4*Hn3/(c3*c4-c3^2)+c3*Hn4/(c4*c3-c4^2))/dt,(-2*Hn3/(c4*c3-c3^2)-2*Hn4/(c4*c3-c4^2))/dt^2],tol,1,2);
        sc=Atol+max(abs(u0),abs(u_new))*Rtol;
        err_coeff = -2*Hn2/(c2*(c2-c4)) ...
                    -2*Hn3/(c3*(c4-c3)) ...
                    +(-2/(c4*(c3-c4)) - 2/(c4*(c4-c2)))*Hn4;
        Er=phipm_simul_iom(dt,Jn,[zero,zero,zero,zero,err_coeff/dt^2],tol,1,2);
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
expHR54s4a_sol=u0;
end

     
% function [time, TDexpR54s4_sol, h, CPU, acc, rej] = TDexpR54s4(F,A,J,t0,t_end,u0,Atol,Rtol,tol,NTS)
% % TDexpR54s4
% % Adaptive 5th-order two-derivative exponential Rosenbrock method
% % with embedded 4th-order estimator.
% %
% % Inputs:
% %   F      - nonlinear/right-hand-side function handle
% %   A      - constant linear part
% %   J      - Jacobian of the nonlinear part
% %   t0     - initial time
% %   t_end  - final time
% %   u0     - initial value
% %   Atol   - absolute tolerance for time stepping
% %   Rtol   - relative tolerance for time stepping
% %   tol    - tolerance passed to phipm_simul_iom
% %   NTS    - initial number of time steps
% %
% % Outputs:
% %   time             - accepted time points
% %   TDexpR54s4_sol   - final solution
% %   h                - accepted step sizes
% %   CPU              - CPU time
% %   acc              - number of accepted steps
% %   rej              - number of rejected steps
% 
%     m    = length(u0);
%     zero = zeros(m,1);
% 
%     % Initial step size
%     dt  = (t_end - t0)/NTS;
%     pow = 1/5;
% 
%     % Nodes
%     c2 = 0.5;
%     c3 = 0.3;
%     c4 = (5*c3 - 3)/(10*c3 - 5);   
% 
%     % Initialization
%     time = t0;
%     h    = [];
%     acc  = 0;
%     rej  = 0;
%     t    = t0;
%     i    = 1;
% 
%     tic
%     while t < t_end
% 
%         % Prevent overshooting final time
%         if t + dt > t_end
%             dt = t_end - t;
%         end
% 
%         % Evaluate at current step
%         Fu0 = F(u0);
%         Jn  =@(u) A*u + J(u0,u);
% 
%         % ----- Stage U_n2 -----
%         % U_n2 = u_n + c2*h*phi_1(c2*h*J_n) F(u_n)
%         Un2 = u0 + phipm_simul_iom(c2*dt, Jn, [zero, Fu0], tol, 1, 2);
% 
%         FUn2 = F(Un2);
%         Hn2  = (J(Un2,FUn2) - J(u0,FUn2));
% 
%         % ----- Stages U_n3 and U_n4 -----
%         % U_ni = u_n + c_i*h*phi_1(c_i*h*J_n)F(u_n)
%         %      + (c_i^3/c2) h^2 phi_3(c_i*h*J_n) H_n2,  i=3,4
%         Utmp = phipm_simul_iom([c3,c4]*dt, Jn, ...
%                [zero, Fu0, zero, (1/c2)*Hn2/dt], tol, 1, 2);
% 
%         Un3 = u0 + Utmp(:,1);
%         Un4 = u0 + Utmp(:,2);
% 
%         FUn3 = F(Un3);
%         FUn4 = F(Un4);
% 
%         Hn3 = (J(Un3,FUn3) - J(u0,FUn3));
%         Hn4 = (J(Un4,FUn4) - J(u0,FUn4));
%         % c = [c3,c4];
%         % Un = zeros(m,2);
%         % Hn = zeros(m,2);
%         % maxNumCompThreads(1);
%         % parfor j=1:2
%         %     Un(:,j) = u0+phipm_simul_iom(c(j)*dt,Jn,[zero,Fu0,zero,1/c2*Hn2/dt],tol,1,2);
%         %     F1 = F(Un(:,j));
%         %     Hn(:,j)=(J(Un(:,j),F1)-J(u0,F1));
%         % end
%         % Hn3 = Hn(:,1);
%         % Hn4 = Hn(:,2);
% 
%         % ----- Main 5th-order update -----
%         coeff_phi3 = c4*Hn3/(c3*(c4-c3)) + c3*Hn4/(c4*(c3-c4));
%         coeff_phi4 = -2*Hn3/(c3*(c4-c3)) - 2*Hn4/(c4*(c3-c4));
% 
%         u_new = u0 + phipm_simul_iom(dt, Jn, ...
%                 [zero, Fu0, zero, coeff_phi3/dt, coeff_phi4/dt^2], tol, 1, 2);
% 
%         % ----- Embedded error estimator -----
%         % E_r = u_{n+1} - \bar{u}_{n+1}
%         err_coeff = -2*Hn2/(c2*(c2-c4)) ...
%                     -2*Hn3/(c3*(c4-c3)) ...
%                     +(-2/(c4*(c3-c4)) - 2/(c4*(c4-c2)))*Hn4;
% 
%         Er = phipm_simul_iom(1, Jn, ...
%              [zero, zero, zero, zero, err_coeff/dt^2], tol, 1, 2);
% 
%         % Scaling for adaptive control
%         sc    = Atol + max(abs(u0), abs(u_new))*Rtol;
%         error = norm(Er./sc, 2)/sqrt(m);
% 
%         % ----- Accept / reject -----
%         if error <= 1
%             t         = t + dt;
%             time(i+1) = t;
%             h(i)      = dt;
%             u0        = u_new;
%             i         = i + 1;
%             acc       = acc + 1;
%         else
%             rej = rej + 1;
%         end
% 
%         % ----- Step size update -----
%         if error < 1e-16
%             s = 2.0;
%         else
%             s = 0.8 * (1/error)^pow;
%         end
% 
%         % Limit step size change
%         dt = dt * min(1.5, max(0.5, s));
%     end
% 
%     CPU = toc;
%     TDexpR54s4_sol = u0;
% end
% function [time,TDexpR54s4_sol,h,CPU,acc,rej] = TDexpR54s4a(F,A,J,t0,t_end,u0,Atol,Rtol,tol,NTS)
% % TDexpR54s4a
% % Adaptive 5th-order two-derivative exponential Rosenbrock method
% % with embedded 4th-order estimator and relaxed stage replacement
% %
% %   (c_i^3/c2) h^2 phi_3(c_i h J_n) H_n2
% % replaced by
% %   (c_i^3/(6*c2)) h^2 H_n2,   i=3,4.
% %
% % Inputs:
% %   F      - nonlinear/right-hand-side function handle
% %   A      - constant linear part
% %   J      - Jacobian of the nonlinear part
% %   t0     - initial time
% %   t_end  - final time
% %   u0     - initial value
% %   Atol   - absolute tolerance for time stepping
% %   Rtol   - relative tolerance for time stepping
% %   tol    - tolerance passed to phipm_simul_iom
% %   NTS    - initial number of time steps
% %
% % Outputs:
% %   time             - accepted time points
% %   TDexpR54s4_sol   - final solution
% %   h                - accepted step sizes
% %   CPU              - CPU time
% %   acc              - number of accepted steps
% %   rej              - number of rejected steps
% 
%     m    = length(u0);
%     zero = zeros(m,1);
% 
%     % Initial step size
%     dt  = (t_end - t0)/NTS;
%     pow = 1/5;
% 
%     % Nodes
%     c2 = 0.3;
%     c3 = 0.3;
%     c4 = (5*c3 - 3)/(10*c3 - 5);
% 
%     % Initialization
%     time = t0;
%     h    = [];
%     acc  = 0;
%     rej  = 0;
%     t    = t0;
%     i    = 1;
% 
%     tic
%     while t < t_end
% 
%         % Prevent overshooting final time
%         if t + dt > t_end
%             dt = t_end - t;
%         end
% 
%         % Current evaluations
%         Ju0 = J(u0);
%         Fu0 = F(u0);
%         Jn  = A + Ju0;
% 
%         % ----- U_n2 and phi1 parts for U_n3, U_n4 -----
%         U = phipm_simul_iom([c2,c4]*dt, Jn, [zero,Fu0], tol, 1, 2);
% 
%         Un2  = u0 + U(:,1);
%         FUn2 = F(Un2);
%         Hn2  = (J(Un2) - Ju0) * FUn2;
% 
%         % ----- Relaxed replacements in U_n3 and U_n4 -----
%         Un3  = u0 + U(:,1) + dt^2*(c3^3/(6*c2))*Hn2;
%         FUn3 = F(Un3);
%         Hn3  = (J(Un3) - Ju0) * FUn3;
% 
%         Un4  = u0 + U(:,2) + dt^2*(c4^3/(6*c2))*Hn2;
%         FUn4 = F(Un4);
%         Hn4  = (J(Un4) - Ju0) * FUn4;
% 
%         % ----- Main 5th-order update -----
%         coeff_phi3 = c4*Hn3/(c3*(c4-c3)) + c3*Hn4/(c4*(c3-c4));
%         coeff_phi4 = -2*Hn3/(c3*(c4-c3)) - 2*Hn4/(c4*(c3-c4));
% 
%         u_new = u0 + dt^2 * phipm_simul_iom(1, dt*Jn, ...
%                 [zero, Fu0/dt, zero, coeff_phi3, coeff_phi4], tol, 1, 2);
% 
%         % ----- Embedded error estimator -----
%         % E_r = u_{n+1} - \bar{u}_{n+1}
%         err_coeff = -2*Hn2/(c2*(c2-c4)) ...
%                     -2*Hn3/(c3*(c4-c3)) ...
%                     +(-2/(c4*(c3-c4)) - 2/(c4*(c4-c2)))*Hn4;
% 
%         Er = dt^2 * phipm_simul_iom(1, dt*Jn, ...
%              [zero, zero, zero, zero, err_coeff], tol, 1, 2);
% 
%         % Error scaling
%         sc    = Atol + max(abs(u0), abs(u_new))*Rtol;
%         error = norm(Er./sc, 2)/sqrt(m);
% 
%         % ----- Accept / reject -----
%         if error <= 1
%             t         = t + dt;
%             time(i+1) = t;
%             h(i)      = dt;
%             u0        = u_new;
%             i         = i + 1;
%             acc       = acc + 1;
%         else
%             rej = rej + 1;
%         end
% 
%         % ----- Step size update -----
%         if error < 1e-16
%             s = 2.0;
%         else
%             s = 0.8 * (1/error)^pow;
%         end
% 
%         % Limit step size change
%         dt = dt * min(1.5, max(0.5, s));
%     end
% 
%     CPU = toc;
%     TDexpR54s4_sol = u0;
% end