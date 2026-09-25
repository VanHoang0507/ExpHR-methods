function [time,expHR54s5_sol,h,CPU,acc,rej] = expHR54s5(F,A,J,t0,t_end,u0,Atol,Rtol,tol,NTS)
% TDexpR54s5
% Adaptive 5th-order two-derivative exponential Rosenbrock method
% with embedded 4th-order estimator.
%
% Inputs:
%   F      - nonlinear/right-hand-side function handle
%   A      - constant linear part
%   J      - Jacobian of the nonlinear part
%   t0     - initial time
%   t_end  - final time
%   u0     - initial value
%   Atol   - absolute tolerance for time stepping
%   Rtol   - relative tolerance for time stepping
%   tol    - tolerance passed to phipm_simul_iom
%   NTS    - initial number of time steps
%
% Outputs:
%   time             - accepted time points
%   TDexpR54s5_sol   - final solution
%   h                - accepted step sizes
%   CPU              - CPU time
%   acc              - number of accepted steps
%   rej              - number of rejected steps

    m    = length(u0);
    zero = zeros(m,1);

    % Initial step size
    dt  = (t_end - t0)/NTS;
    pow = 1/5;

    % Nodes
    c2 = 0.3;
    c3 = c2;
    c4 = 0.9;
    c5 = 1.0;

    % Coefficients rho_{ji}
    % For i in {3,4,5}, with the other two nodes denoted by k,l:
    % rho_{3i} = c_k*c_l / ( c_i (c_i-c_l)(c_i-c_k) )
    % rho_{4i} = -2(c_k+c_l) / ( c_i (c_i-c_l)(c_i-c_k) )
    % rho_{5i} = 6 / ( c_i (c_i-c_l)(c_i-c_k) )

    rho33 =  c4*c5 / (c3*(c3-c4)*(c3-c5));
    rho34 = -2*(c4+c5) / (c3*(c3-c4)*(c3-c5));
    rho35 =  6 / (c3*(c3-c4)*(c3-c5));

    rho43 =  c3*c5 / (c4*(c4-c3)*(c4-c5));
    rho44 = -2*(c3+c5) / (c4*(c4-c3)*(c4-c5));
    rho45 =  6 / (c4*(c4-c3)*(c4-c5));

    rho53 =  c3*c4 / (c5*(c5-c3)*(c5-c4));
    rho54 = -2*(c3+c4) / (c5*(c5-c3)*(c5-c4));
    rho55 =  6 / (c5*(c5-c3)*(c5-c4));

    % Initialization
    time = t0;
    h    = [];
    acc  = 0;
    rej  = 0;
    t    = t0;
    i    = 1;

    tic
    while t < t_end

        % Prevent overshooting final time
        if t + dt > t_end
            dt = t_end - t;
        end

        % Current evaluations
        Fu0 = F(u0);
        Jn  = @(u) A*u + J(u0,u);

        % ----- Stage U_n2 -----
        % U_n2 = u_n + c2*h*phi_1(c2*h*J_n)F(u_n)
        Un2  = u0 + phipm_simul_iom(c2*dt, Jn, [zero, Fu0], tol, 1, 2);
        FUn2 = F(Un2);
        Hn2  = (J(Un2,FUn2) - J(u0,FUn2));

        % ----- Stages U_n3, U_n4, U_n5 -----
        % U_ni = u_n + c_i*h*phi_1(c_i*h*J_n)F(u_n)
        %      + (c_i^3/c2) h^2 phi_3(c_i*h*J_n)H_n2,  i=3,4,5

        U1=phipm_simul_iom([c3,c4,c5]*dt,Jn,[zero,Fu0,zero,1/c2*Hn2/dt],tol,1,2);
        Un3=u0+U1(:,1);
        FUn3 = F(Un3);
        Hn3=(J(Un3,FUn3)-J(u0,FUn3));
        Un4=u0+U1(:,2);
        FUn4 = F(Un4);
        Hn4=(J(Un4,FUn4)-J(u0,FUn4));
        Un5=u0+U1(:,3);
        FUn5 = F(Un5); 
        Hn5=(J(Un5,FUn5)-J(u0,FUn5));
        
% c = [c3,c4,c5];
% Hn = zeros(length(u0),3);
% maxNumCompThreads(1);
% parfor j = 1:3
%     Unj  = u0 + phipm_simul_iom(c(j)*dt, Jn, [zero,Fu0,zero,1/c2*Hn2/dt], tol, 1, 2);
%     FUnj = F(Unj);
%     Hn(:,j) = J(Unj,FUnj) - J(u0,FUnj);
% end
% 
% Hn3 = Hn(:,1);
% Hn4 = Hn(:,2);
% Hn5 = Hn(:,3);

        % ----- Main 5th-order update -----
        coeff_phi3 = rho33*Hn3 + rho43*Hn4 + rho53*Hn5;
        coeff_phi4 = rho34*Hn3 + rho44*Hn4 + rho54*Hn5;
        coeff_phi5 = rho35*Hn3 + rho45*Hn4 + rho55*Hn5;

        u_new = u0 + phipm_simul_iom(dt, Jn, ...
                [zero, Fu0, zero, coeff_phi3/dt, coeff_phi4/dt^2, coeff_phi5/dt^3], tol, 1, 2);

        % ----- Embedded 4th-order estimator -----
        % E_r = u_{n+1} - \bar{u}_{n+1} = h^2 phi_5(hJ_n) * sum_{i=3}^5 rho_{5i} H_{ni}
        Er =  phipm_simul_iom(dt, Jn, ...
             [zero, zero, zero, zero, zero, coeff_phi5/dt^3], tol, 1, 2);

        % Error scaling
        sc    = Atol + max(abs(u0), abs(u_new))*Rtol;
        error = norm(Er./sc, 2)/sqrt(m);

        % ----- Accept / reject -----
        if error <= 1
            t         = t + dt;
            time(i+1) = t;
            h(i)      = dt;
            u0        = u_new;
            i         = i + 1;
            acc       = acc + 1;
        else
            rej = rej + 1;
        end

        % ----- Step size update -----
        if error < 1e-16
            s = 2.0;
        else
            s = 0.8 * (1/error)^pow;
        end

        % Limit step size change
        dt = dt * min(1.5, max(0.5, s));
    end

    CPU = toc;
    expHR54s5_sol = u0;
end