clc
clear all
close all
addpath('../integrators','../phipmsimuliom','./Data')
outFolder = 'Test Results';
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end
global a b nx ny kappa;
% Parameters

nx  = 100;             % number of grid points in x
ny  = 100;             % number of grid points in y

Du = 0.05;             % D1
Dv = 1.0;              % D2

a  = 0.1305;
b  = 0.7695;
kappa = 100;

Lx = 1;
Ly = 1;

dx = Lx/(nx-1);
dy = Ly/(ny-1);

tspan = [0, 1];       % Time range

% Initial condition
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
Dx2(1,1)   = -2;  Dx2(1,2)   = 2;
Dx2(nx,nx) = -2;  Dx2(nx,nx-1) = 2;
Dx2 = Dx2/dx^2;

ey = ones(ny,1);
Dy2 = spdiags([ey -2*ey ey],[-1 0 1],ny,ny);
Dy2(1,1)   = -2;  Dy2(1,2)   = 2;
Dy2(ny,ny) = -2;  Dy2(ny,ny-1) = 2;
Dy2 = Dy2/dy^2;

%% 2D Laplacian
L2 = kron(speye(ny),Dx2) + kron(Dy2,speye(nx));

%% Block linear operator
Iu = [1,0;0,0];
Iv = [0,0;0,1];

A = kron(Iu,Du*L2) + kron(Iv,Dv*L2);

g  = @(Y) rhs(Y,nx,ny,a,b,kappa);
F  = @(Y) A*Y + g(Y);
Jg = @(Y,X) Jac_g_vec(Y,X,nx,ny,kappa);

load ref_sol_Schnakenberg_k100_t1.mat ode15s_sol_end
exac_sol = ode15s_sol_end;
% %This part is the comparison between embedded schemes
 
 tol_expHR42s2 = [10^-2.46 10^-4 10^-5 10^-6.24 10^-9.2]; % N = 150
 tol_expHR43s3 = [10^-3.05 10^-4.1 10^-4.88 10^-6.18]; 
 tol_expRB43s3 = [10^-2.73 10^-3.65 10^-4.6 10^-5.5];
 tol_pexpRB43s3 = [10^-3 10^-3.95 10^-4.99 10^-5.9];
 tol_expRB42s2 = [10^-2.51 10^-3.22 10^-4.55 10^-6.3];
 
% % % % % 
for i=1:4

    [~,expHR42s2_sol,~,CPU,~,~]=expHR32_p(F,A,Jg,tspan(1),tspan(2),Y0,tol_expHR42s2(i),tol_expHR42s2(i),1e-10,1/2,150);
    expHR42s2_err(i)=norm(exac_sol-expHR42s2_sol,'inf');
    CPU_expHR42s2(i)=CPU;

    [~,expHR43s3_sol,~,CPU,~,~]=expHR43s3(F,A,Jg,tspan(1),tspan(2),Y0,tol_expHR43s3(i),tol_expHR43s3(i),1e-10,150);
    expHR43s3_err(i)=norm(exac_sol-expHR43s3_sol,'inf');
    CPU_TDexpR43s3(i)=CPU;

    [~,expRB42s2_sol,~,CPU,~,~]=expRB42s2(F,A,g,Jg,tspan(1),tspan(2),Y0,tol_expRB42s2(i),tol_expRB42s2(i),1e-10,150);
    expRB42s2_err(i)=norm(exac_sol-expRB42s2_sol,'inf');
    CPU_expRB42s2(i)=CPU;

    [~,pexpRB43s3_sol,~,CPU,~,~]=pexpRB43s3(F,A,g,Jg,tspan(1),tspan(2),Y0,tol_pexpRB43s3(i),tol_pexpRB43s3(i),1e-10,150);  
    pexpRB43s3_err(i)=norm(exac_sol-pexpRB43s3_sol,'inf');
    CPU_pexpRB43s3(i)=CPU;

    [~,expRB43s3_sol,~,CPU,~,~]=expRB43s3(F,A,g,Jg,tspan(1),tspan(2),Y0,tol_expRB43s3(i),tol_expRB43s3(i),1e-10,150);   
    expRB43s3_err(i)=norm(exac_sol-expRB43s3_sol,'inf');
    CPU_expRB43s3(i)=CPU;

end

figure(1)
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 6.8 5.2]);
set(gcf,'PaperPositionMode','auto');

semilogy(CPU_expHR42s2,expHR42s2_err,'o-','LineWidth',1, ...
    'Color',[0.8660 0.3290 0.0000],'MarkerSize',10); hold on      % TDexpRB42s2 orange circle

semilogy(CPU_TDexpR43s3,expHR43s3_err,'v-','LineWidth',1, ...
    'Color',[0.9290 0.6940 0.1250],'MarkerSize',10);              % TDexpRB43s3 yellow down-triangle

semilogy(CPU_expRB43s3,expRB43s3_err,'d-','LineWidth',1, ...
    'Color',[0.8660 0.3290 0.0000],'MarkerSize',10);              % exprb43 orange diamond

semilogy(CPU_pexpRB43s3,pexpRB43s3_err,'s-','LineWidth',1, ...
    'Color',[0.9290 0.6940 0.1250],'MarkerSize',10);              % pexprb43 yellow square

semilogy(CPU_expRB42s2,expRB42s2_err,'^-','LineWidth',1, ...
    'Color',[0.0660 0.4430 0.7450],'MarkerSize',10);              % exprb42 blue up-triangle

xlabel('CPU time','FontSize',10)
ylabel('Error','FontSize',10)
grid on
box on

ax = gca;
ax.Units = 'normalized';
ax.Position = [0.24 0.16 0.72 0.72];
ax.FontSize = 10;
set(ax,'TickLength',2*get(ax,'TickLength'))

set(ax,'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2])

axis square

lgd = legend('expHR42s2','expHR43s3','exprb43','pexprb43','exprb42s2', ...
    'Location','westoutside');
lgd.FontSize = 10;


baseName = 'expHR_Schnakenberg_CPU4_embedded';
print(gcf,fullfile(outFolder, baseName), '-depsc');
saveas(gcf,fullfile(outFolder, [baseName '.png']));

rmpath('../integrators','../phipmsimuliom','./Data')

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