clc
clear all
close all

addpath('../integrators','../phipmsimuliom')

global a b nx ny

%% Parameters
nx  = 100;             % number of grid points in x
ny  = 100;             % number of grid points in y

Du = 0.05;             % D1
Dv = 1.0;              % D2

a  = 0.1305;
b  = 0.7695;

kappa_list = [100, 1000];   % values of kappa to test

Lx = 1;
Ly = 1;

dx = Lx/(nx-1);
dy = Ly/(ny-1);

tspan = [0, 1];        % Time range

%% Create output folder
outFolder = 'Test Results';
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

%% Initial condition
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[X,Y] = meshgrid(x,y);
X = X';
Y = Y';

initial_condition_u = a + b + 1e-3*exp(-100*((X-1/3).^2 + (Y-1/2).^2));
initial_condition_v = b/(a+b)^2 * ones(nx,ny);

U0 = reshape(initial_condition_u,nx*ny,1);
V0 = reshape(initial_condition_v,nx*ny,1);
Y0 = [U0;V0];

%% 1D Neumann Laplacians
ex = ones(nx,1);
Dx2 = spdiags([ex -2*ex ex],[-1 0 1],nx,nx);
Dx2(1,1)   = -2;  
Dx2(1,2)   = 2;
Dx2(nx,nx) = -2;  
Dx2(nx,nx-1) = 2;
Dx2 = Dx2/dx^2;

ey = ones(ny,1);
Dy2 = spdiags([ey -2*ey ey],[-1 0 1],ny,ny);
Dy2(1,1)   = -2;  
Dy2(1,2)   = 2;
Dy2(ny,ny) = -2;  
Dy2(ny,ny-1) = 2;
Dy2 = Dy2/dy^2;

%% 2D Laplacian
L2 = kron(speye(ny),Dx2) + kron(Dy2,speye(nx));

%% Block linear operator
Iu = [1,0;0,0];
Iv = [0,0;0,1];
A = kron(Iu,Du*L2) + kron(Iv,Dv*L2);

%% Target output times
tplot = [0, 0.5, 1];

%% Loop over kappa values
for kk = 1:length(kappa_list)

    kappa = kappa_list(kk);

    g  = @(Y) rhs(Y,nx,ny,a,b,kappa);
    F  = @(Y) A*Y + g(Y);
    Jg = @(Y,X) Jac_g_vec(Y,X,nx,ny,kappa);

    % Solve and store snapshots
    [time,Sol,h,CPU,acc,rej] = expHR54s4_plot(F,A,Jg,tspan(1),tspan(2),Y0,10^-7.5,10^-7.5,1e-10,100);

    fprintf('kappa = %d: CPU = %.4f, accepted = %d, rejected = %d\n', ...
        kappa, CPU, acc, rej);

    % Extract snapshots nearest to desired times
    idx0  = find_closest_time_index(time,0);
    idx05 = find_closest_time_index(time,0.5);
    idx1  = find_closest_time_index(time,1);

    X0  = Sol(:,idx0);
    X05 = Sol(:,idx05);
    X1  = Sol(:,idx1);

    % ----- u component -----
    u0_plot  = reshape(X0(1:nx*ny), [nx ny]);
    u05_plot = reshape(X05(1:nx*ny), [nx ny]);
    u1_plot  = reshape(X1(1:nx*ny), [nx ny]);

    fig = figure('Position',[100 100 1800 550]);

    ax1 = subplot(1,3,1);
    surf(x,y,u0_plot')
    shading interp
    view(2)
    axis([0 1 0 1])
    axis square
    title('u at t = 0')
    xlabel('x')
    ylabel('y')
    colorbar

    ax2 = subplot(1,3,2);
    surf(x,y,u05_plot')
    shading interp
    view(2)
    axis([0 1 0 1])
    axis square
    title('u at t = 0.5')
    xlabel('x')
    ylabel('y')
    colorbar

    ax3 = subplot(1,3,3);
    surf(x,y,u1_plot')
    shading interp
    view(2)
    axis([0 1 0 1])
    axis square
    title('u at t = 1')
    xlabel('x')
    ylabel('y')
    colorbar

    sgtitle(sprintf('\\kappa = %d',kappa),'FontSize',16,'FontWeight','bold')

    set(ax1,'Position',[0.05 0.16 0.23 0.72])
    set(ax2,'Position',[0.37 0.16 0.23 0.72])
    set(ax3,'Position',[0.69 0.16 0.23 0.72])

    % Save figure
    baseName = sprintf('expHR_Schnakenberg_kappa_%d_plot',kappa);
    print(fig, fullfile(outFolder, baseName), '-depsc');
    saveas(fig, fullfile(outFolder, [baseName '.png']));

end

%% =======================================================================
function idx = find_closest_time_index(time,tt)
    [~,idx] = min(abs(time-tt));
end

function g = rhs(Y,nx,ny,a,b,kappa)

    N = nx*ny;

    U = Y(1:N);
    V = Y(N+1:end);

    U2V = (U.^2).*V;

    gu = kappa*(a - U + U2V);
    gv = kappa*(b - U2V);

    g = [gu; gv];
end

function JX = Jac_g_vec(Y, X, nx, ny, kappa)

    N = nx*ny;

    U = Y(1:N);
    V = Y(N+1:end);

    UV = U.*V;
    U2 = U.^2;

    JX = [kappa*(-1 + 2*UV).*X(1:N) + kappa*U2.*X(N+1:2*N);
         -2*kappa*UV.*X(1:N) - kappa*U2.*X(N+1:2*N)];
end

function [time,expHR54s4_sol,h,CPU,acc,rej] = expHR54s4_plot(F,A,J,t0,t_end,u0,Atol,Rtol,tol,NTS)

m = length(u0);
zero = zeros(m,1);
dt = (t_end-t0)/NTS;
pow = 1/5;
c2 = 0.1;
c3 = 0.7;
c4 = (c3*5-3)/(10*c3-5);

time = t0;
i = 1;
expHR54s4_sol(:,1) = u0;
h(1) = dt;
acc = 0;
rej = 0;
t = t0;

tic
while t < t_end

    if t + dt > t_end
        dt = t_end - t;
    end

    Fu0 = F(u0);
    Jn = @(u) A*u + J(u0,u);

    Un2 = u0 + phipm_simul_iom(c2*dt,Jn,[zero,Fu0],tol,1,2);
    FUn2 = F(Un2);
    Hn2 = (J(Un2,FUn2) - J(u0,FUn2));

    U1 = phipm_simul_iom([c4,c3]*dt,Jn,[zero,Fu0,zero,(1/c2)*Hn2/dt],tol,1,2);

    Un3 = u0 + U1(:,2);
    FUn3 = F(Un3);
    Hn3 = (J(Un3,FUn3) - J(u0,FUn3));

    Un4 = u0 + U1(:,1);
    FUn4 = F(Un4);
    Hn4 = (J(Un4,FUn4) - J(u0,FUn4));

    u_new = u0 + phipm_simul_iom(dt,Jn,...
        [zero,Fu0,zero,...
        (c4*Hn3/(c3*c4-c3^2)+c3*Hn4/(c4*c3-c4^2))/dt,...
        (-2*Hn3/(c4*c3-c3^2)-2*Hn4/(c4*c3-c4^2))/dt^2],tol,1,2);

    sc = Atol + max(abs(u0),abs(u_new))*Rtol;

    err_coeff = -2*Hn2/(c2*(c2-c4)) ...
              - 2*Hn3/(c3*(c4-c3)) ...
              + (-2/(c4*(c3-c4)) - 2/(c4*(c4-c2)))*Hn4;

    Er = phipm_simul_iom(dt,Jn,[zero,zero,zero,zero,err_coeff/dt^2],tol,1,2);

    error = norm(Er./sc,2)/sqrt(m);

    if error <= 1
        t = t + dt;
        i = i + 1;
        time(i) = t;
        h(i) = dt;
        u0 = u_new;
        expHR54s4_sol(:,i) = u0;
        acc = acc + 1;
    else
        rej = rej + 1;
    end

    if error < 1e-16
        s = 2;
    else
        s = 0.8 * (1/error)^pow;
    end

    dt = dt * min(1.5, max(0.5, s));
end

CPU = toc;

end