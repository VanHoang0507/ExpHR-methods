clc
clear all
close all

addpath('../integrators','../phipmsimuliom','./Data')

%% ========================================================================
% Parameters
% ========================================================================

nx  = 100;
ny  = 100;

Du = 0.05;
Dv = 1.0;

a = 0.1305;
b = 0.7695;

% Test both stiffness parameters
kappa_list = [100, 1000];

Lx = 1;
Ly = 1;

dx = Lx/(nx-1);
dy = Ly/(ny-1);

tspan = [0,1];

%% ========================================================================
% Output folder
% ========================================================================

output_folder = 'Test Results';

if ~exist(output_folder,'dir')
    mkdir(output_folder);
end

%% ========================================================================
% Initial condition
% ========================================================================

x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);

[X,Y] = meshgrid(x,y);

X = X';
Y = Y';

initial_condition_u = ...
    a + b + 1e-3*exp(-100*((X-1/3).^2 + (Y-1/2).^2));

initial_condition_v = ...
    b/(a+b)^2 * ones(nx,ny);

U0 = reshape(initial_condition_u,nx*ny,1);
V0 = reshape(initial_condition_v,nx*ny,1);

Y0 = [U0;V0];

%% ========================================================================
% 1D Neumann Laplacians
% ========================================================================

ex = ones(nx,1);

Dx2 = spdiags([ex -2*ex ex],[-1 0 1],nx,nx);

Dx2(1,1)     = -2;
Dx2(1,2)     =  2;

Dx2(nx,nx)   = -2;
Dx2(nx,nx-1) =  2;

Dx2 = Dx2/dx^2;


ey = ones(ny,1);

Dy2 = spdiags([ey -2*ey ey],[-1 0 1],ny,ny);

Dy2(1,1)     = -2;
Dy2(1,2)     =  2;

Dy2(ny,ny)   = -2;
Dy2(ny,ny-1) =  2;

Dy2 = Dy2/dy^2;

%% ========================================================================
% 2D Laplacian
% ========================================================================

L2 = kron(speye(ny),Dx2) + ...
     kron(Dy2,speye(nx));

%% ========================================================================
% Block linear operator
% ========================================================================

Iu = [1,0;
      0,0];

Iv = [0,0;
      0,1];

A = kron(Iu,Du*L2) + ...
    kron(Iv,Dv*L2);

%% ========================================================================
% Storage for results
% ========================================================================

Results = struct();

%% ========================================================================
% LOOP OVER KAPPA
% ========================================================================

for ik = 1:length(kappa_list)

    kappa = kappa_list(ik);

    fprintf('\n');
    fprintf('=============================================\n');
    fprintf('Running Schnakenberg problem: kappa = %d\n',kappa);
    fprintf('=============================================\n');

    %% --------------------------------------------------------------------
    % Nonlinear function and Jacobian
    % ---------------------------------------------------------------------

    g   = @(Y) rhs(Y,nx,ny,a,b,kappa);

    F   = @(Y) A*Y + g(Y);

    Jg  = @(Y,X) ...
        Jac_g_vec(Y,X,nx,ny,kappa);

    Jg1 = @(Y) ...
        Jac_g(Y,nx,ny,kappa);

    %% --------------------------------------------------------------------
    % Select reference solution and tolerances
    % ---------------------------------------------------------------------

    switch kappa

        case 100

            ref_file = ...
                'ref_sol_Schnakenberg_k100_t1.mat';

            tol_expHR54s4 = 10^-3.2009;

            tol_ode15s = 10^-5.69;


        case 1000

            ref_file = ...
                'ref_sol_Schnakenberg_k1000_t1.mat';

            % More stringent tolerances for the stiffer problem
            tol_expHR54s4 = 10^-10;
            % previous possible value: 10^-8.74

            tol_ode15s = 10^-12.899;

    end

    fprintf('Reference file : %s\n',ref_file);
    fprintf('expHR54s4 tol  : %.4e\n',tol_expHR54s4);
    fprintf('ode15s tol     : %.4e\n\n',tol_ode15s);

    %% --------------------------------------------------------------------
    % Load reference solution
    % ---------------------------------------------------------------------

    ref_data = load(ref_file,'ode15s_sol_end');

    exac_sol = ref_data.ode15s_sol_end;

    %% --------------------------------------------------------------------
    % ode15s options
    % ---------------------------------------------------------------------

    opts = odeset( ...
        'RelTol',tol_ode15s, ...
        'AbsTol',tol_ode15s, ...
        'Jacobian',@(t,Y) A + Jg1(Y), ...
        'InitialStep',0.01, ...
        'Stats','on');

    %% ====================================================================
    % expHR54s4
    % =====================================================================

    fprintf('Running expHR54s4...\n');

    [t_expHR54s4, ...
     expHR54s4_sol, ...
     h_expHR54s4, ...
     CPU_expHR54s4, ...
     acc_expHR54s4, ...
     rej_expHR54s4] = ...
        expHR54s4( ...
        F, ...
        A, ...
        Jg, ...
        tspan(1), ...
        tspan(2), ...
        Y0, ...
        tol_expHR54s4, ...
        tol_expHR54s4, ...
        1e-10, ...
        100);

    %% Error

    err_expHR54s4 = ...
        norm(exac_sol - expHR54s4_sol,'inf');

    fprintf('\n');
    fprintf('expHR54s4 results\n');
    fprintf('------------------------------\n');
    fprintf('Error     = %.6e\n',err_expHR54s4);
    fprintf('Accepted  = %d\n',acc_expHR54s4);
    fprintf('Rejected  = %d\n',rej_expHR54s4);
    fprintf('CPU       = %.6f sec\n',CPU_expHR54s4);

    %% ====================================================================
    % ode15s
    % =====================================================================

    fprintf('\nRunning ode15s...\n');

    tic

    [t_ode15s,ode15s_sol] = ...
        ode15s( ...
        @(t,Y) A*Y + g(Y), ...
        tspan, ...
        Y0, ...
        opts);

    CPU_ode15s = toc;

    ode15s_sol = ode15s_sol';

    ode15s_sol_end = ode15s_sol(:,end);

    %% Error

    err_ode15s = ...
        norm(exac_sol - ode15s_sol_end,'inf');

    fprintf('\n');
    fprintf('ode15s results\n');
    fprintf('------------------------------\n');
    fprintf('Error = %.6e\n',err_ode15s);
    fprintf('CPU   = %.6f sec\n',CPU_ode15s);

    %% ====================================================================
    % Store results
    % =====================================================================

    Results(ik).kappa = kappa;

    Results(ik).tol_expHR54s4 = tol_expHR54s4;
    Results(ik).tol_ode15s = tol_ode15s;

    Results(ik).err_expHR54s4 = err_expHR54s4;
    Results(ik).err_ode15s = err_ode15s;

    Results(ik).CPU_expHR54s4 = CPU_expHR54s4;
    Results(ik).CPU_ode15s = CPU_ode15s;

    Results(ik).acc_expHR54s4 = acc_expHR54s4;
    Results(ik).rej_expHR54s4 = rej_expHR54s4;

    Results(ik).t_expHR54s4 = t_expHR54s4;
    Results(ik).h_expHR54s4 = h_expHR54s4;

    Results(ik).t_ode15s = t_ode15s;

    %% ====================================================================
    % Step sizes
    % =====================================================================

    h_ode15s = diff(t_ode15s);

    %% ====================================================================
    % Plot
    % =====================================================================

    fig = figure;

    set(fig,'Units','inches');
    set(fig,'Position',[1 1 6.8 5.2]);
    set(fig,'PaperPositionMode','auto');

    semilogy( ...
        t_expHR54s4(1:end-1), ...
        h_expHR54s4(2:end), ...
        '-.', ...
        'LineWidth',2, ...
        'Color',[0.8660 0.3290 0.0000]);

    hold on

    semilogy( ...
        t_ode15s(1:end-1), ...
        h_ode15s, ...
        '-', ...
        'LineWidth',2, ...
        'Color',[0.0660 0.4430 0.7450]);

    xlabel('Time','FontSize',12)

    ylabel('Step size','FontSize',12)

    grid on
    box on

    ax = gca;

    ax.Units = 'normalized';
    ax.Position = [0.24 0.16 0.72 0.72];
    ax.FontSize = 10;

    set(ax, ...
        'TickLength', ...
        2*get(ax,'TickLength'))

    axis square

    lgd = legend( ...
        'expHR54s4', ...
        'ode15s', ...
        'Location','westoutside');

    lgd.FontSize = 10;

    title(sprintf('\\kappa = %d',kappa), ...
        'FontSize',12)

    %% ====================================================================
    % Save figure
    % =====================================================================

    eps_name = fullfile( ...
        output_folder, ...
        sprintf('Stepsize_vs_ode15s_kappa_%d.eps',kappa));

    png_name = fullfile( ...
        output_folder, ...
        sprintf('Stepsize_vs_ode15s_kappa_%d.png',kappa));

    print(fig,eps_name,'-depsc');

    exportgraphics(fig,png_name,'Resolution',300);

end

%% ========================================================================
% Display comparison table
% ========================================================================

fprintf('\n\n');
fprintf('=============================================================\n');
fprintf('                FINAL COMPARISON\n');
fprintf('=============================================================\n\n');

for ik = 1:length(Results)

    fprintf('kappa = %d\n',Results(ik).kappa);

    fprintf('   expHR54s4:\n');
    fprintf('      Error = %.6e\n',Results(ik).err_expHR54s4);
    fprintf('      CPU   = %.6f sec\n',Results(ik).CPU_expHR54s4);
    fprintf('      Acc   = %d\n',Results(ik).acc_expHR54s4);
    fprintf('      Rej   = %d\n',Results(ik).rej_expHR54s4);

    fprintf('   ode15s:\n');
    fprintf('      Error = %.6e\n',Results(ik).err_ode15s);
    fprintf('      CPU   = %.6f sec\n',Results(ik).CPU_ode15s);

    fprintf('\n');

end

%% ========================================================================
% Save all numerical results
% ========================================================================

save( ...
    fullfile(output_folder,'Schnakenberg_expHR54s4_vs_ode15s.mat'), ...
    'Results');

%% ========================================================================
% Nonlinear term
% ========================================================================

function g = rhs(Y,nx,ny,a,b,kappa)

    N = nx*ny;

    U = Y(1:N);
    V = Y(N+1:end);

    U2V = (U.^2).*V;

    gu = kappa*(a - U + U2V);
    gv = kappa*(b - U2V);

    g = [gu;gv];

end

%% ========================================================================
% Full Jacobian
% ========================================================================

function J = Jac_g(Y,nx,ny,kappa)

    N = nx*ny;

    U = Y(1:N);
    V = Y(N+1:end);

    UV = U.*V;
    U2 = U.^2;

    I11 = [1,0;
           0,0];

    I12 = [0,1;
           0,0];

    I21 = [0,0;
           1,0];

    I22 = [0,0;
           0,1];

    J = ...
        kron(I11, ...
        spdiags(kappa*(-1 + 2*UV),0,N,N)) ...
        ...
      + kron(I12, ...
        spdiags(kappa*U2,0,N,N)) ...
        ...
      + kron(I21, ...
        spdiags(-2*kappa*UV,0,N,N)) ...
        ...
      + kron(I22, ...
        spdiags(-kappa*U2,0,N,N));

end

%% ========================================================================
% Jacobian-vector product
% ========================================================================

function JX = Jac_g_vec(Y,X,nx,ny,kappa)

    N = nx*ny;

    U = Y(1:N);
    V = Y(N+1:end);

    UV = U.*V;
    U2 = U.^2;

    Xu = X(1:N);
    Xv = X(N+1:2*N);

    JX = [ ...
        kappa*(-1 + 2*UV).*Xu ...
        + kappa*U2.*Xv;
        ...
        -2*kappa*UV.*Xu ...
        - kappa*U2.*Xv ];

end