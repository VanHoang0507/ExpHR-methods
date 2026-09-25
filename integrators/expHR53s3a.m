function [time, expHR53s3_sol, h, CPU, acc, rej] = expHR53s3a(F,A,J,t0,t_end,u0,Atol,Rtol,tol,NTS)
% TDexpR53s3a
% Adaptive 5th-order two-derivative exponential Rosenbrock method
% with embedded 3rd-order estimator and relaxed stage replacement.
%
% Relaxed replacement in U_n3:
%   c2^2*(b2(0)/b3(0))*h^2*phi_3(c2*h*J_n)H_n2
% + (c3^3/c2)*h^2*phi_3(c3*h*J_n)H_n2
% is replaced by
%   ( c2^2*b2(0)/(6*b3(0)) + c3^3/(6*c2) ) * h^2 * H_n2.
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
%   TDexpR53s3_sol   - final solution
%   h                - accepted step sizes
%   CPU              - CPU time
%   acc              - number of accepted steps
%   rej              - number of rejected steps

    m    = length(u0);
    zero = zeros(m,1);

    % Initial step size
    dt  = (t_end - t0)/NTS;
    pow = 1/4;   % embedded order 3 -> controller exponent 1/4

    %c2, c3
    c2 = 0.3;
    c3 = (5*c2 - 3)/(5*(2*c2 - 1));

    % Coefficients at zero
    b20 = (c3/6 - 1/12) / (c2*(c3-c2));
    b30 = (c2/6 - 1/12) / (c3*(c2-c3));

    % Relaxed stage coefficient
    a320 = c2^2*b20/(6*b30) + c3^3/(6*c2);

    % Main update coefficients
    beta3_2 =  c3/(c2*(c3-c2));
    beta3_3 =  c2/(c3*(c2-c3));

    beta4_2 = -2/(c2*(c3-c2));
    beta4_3 = -2/(c3*(c2-c3));

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

        % Evaluate at current step
        
        Fu0 = F(u0);
        Jn  =@(u) A*u + J(u0,u);

        % ----- Stage U_n2 -----
        U = phipm_simul_iom([c2, c3]*dt, Jn, [zero, Fu0], tol, 1, 2);
        % c = [c2,c3];
        % U = zeros(m,2);
        % maxNumCompThreads(1);
        % parfor j=1:2
        %     U(:,j) = u0+phipm_simul_iom(c(j)*dt,Jn,[zero,Fu0],tol,1,2);
        % end
        Un2  = u0 + U(:,1);
        FUn2 = F(Un2);
        Hn2  = (J(Un2,FUn2) - J(u0,FUn2));

        % ----- Relaxed stage U_n3 -----
        % U_n3 = u_n + c3*h*phi_1(c3*h*J_n)F(u_n) + a320*h^2*H_n2
        Un3  = u0 + U(:,2) + dt^2*a320*Hn2;
        FUn3 = F(Un3);
        Hn3  = (J(Un3,FUn3) - J(u0,FUn3));

        % ----- Main 5th-order update -----
        coeff_phi3 = beta3_2*Hn2 + beta3_3*Hn3;
        coeff_phi4 = beta4_2*Hn2 + beta4_3*Hn3;

        u_new = u0 + phipm_simul_iom(dt, Jn, ...
                [zero, Fu0, zero, coeff_phi3/dt, coeff_phi4/dt^2], tol, 1, 2);

        % ----- Embedded error estimator -----
        % E_r = u_{n+1} - \bar{u}_{n+1}
        Er = phipm_simul_iom(dt, Jn, ...
             [zero, zero, zero, zero, coeff_phi4/dt^2], tol, 1, 2);

        % Scaling for adaptive control
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

        dt = dt * min(1.5, max(0.5, s));
    end

    CPU = toc;
    expHR53s3_sol = u0;
end