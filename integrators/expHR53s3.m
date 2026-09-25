function [time, expHR53s3_sol, h, CPU, acc, rej] = expHR53s3(F,A,J,t0,t_end,u0,Atol,Rtol,tol,NTS)
% TDexpR53s3
% Adaptive 5th-order two-derivative exponential Rosenbrock method
% with embedded 3rd-order estimator.
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

    % c2,c3
    c2=0.3; 
    c3=(c2*5-3)/(10*c2-5);

    % Coefficients at zero
    b20 = (c3/6 - 1/12) / (c2*(c3-c2));
    b30 = (c2/6 - 1/12) / (c3*(c2-c3));

    % Initialization
    time = t0;
    h    = [];
    acc  = 0;
    rej  = 0;
    t    = t0;
    i    = 1;
    maxNumCompThreads(1);
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
        Un2 = u0 + phipm_simul_iom(c2*dt, Jn, [zero, Fu0], tol, 1, 2);

        FUn2 = F(Un2);
        Hn2  = (J(Un2,FUn2) - J(u0,FUn2)) ;

        % ----- Stage U_n3 -----
        % U_n3 = u_n + c3*h*phi1(c3*h*J_n)F(u_n)
        %      + c2^2*(b20/b30)*h^2*phi3(c2*h*J_n)H_n2
        %      + (c3^3/c2)*h^2*phi3(c3*h*J_n)H_n2

        part_c3 = phipm_simul_iom(c3*dt, Jn, ...
                  [zero, Fu0, zero, (1/c2)*Hn2/dt], tol, 1, 2);

        part_c2 = phipm_simul_iom(c2*dt, Jn, ...
                  [zero, zero, zero, (b20/(c2*b30))*Hn2/dt], tol, 1, 2);

        Un3 = u0 + part_c3 + part_c2;

        FUn3 = F(Un3);
        Hn3  = (J(Un3,FUn3) - J(u0,FUn3));

        % ----- Main 5th-order update -----
        coeff_phi3 = c3*Hn2/(c2*(c3-c2)) + c2*Hn3/(c3*(c2-c3));
        coeff_phi4 = -2*Hn2/(c2*(c3-c2)) - 2*Hn3/(c3*(c2-c3));

        u_new = u0 + phipm_simul_iom(dt, Jn, ...
                [zero, Fu0, zero, coeff_phi3/dt, coeff_phi4/dt^2], tol, 1, 2);

        % ----- Embedded error estimator -----
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