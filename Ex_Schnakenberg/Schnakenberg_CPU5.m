clc
clear all
close all
addpath('../integrators','../phipmsimuliom','./Data')

%% Create output folder
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This part is the comparison between embedded schemes
tol_expHR53s3=[10^-4.525 10^-5.34 10^-8.65  10^-9.23]; %N = 100; 
tol_expHR54s4=[10^-4.42 10^-5.295 10^-7.5 10^-8.6]; 
tol_expHR53s3a=[10^-5.43 10^-7 10^-8.65 10^-9.23];
tol_expHR54s4a=[10^-6 10^-7.2 10^-8.2 10^-9.15]; 
tol_expHR54s5=[10^-4.3 10^-6.04 10^-8.9 10^-9.49]; 
tol_expRB53s3 =[10^-5.0 10^-6.5 10^-8.9 10^-10.1]; 
tol_pexpRB54s4=[10^-5.07 10^-6.375 10^-8.58 10^-9.65 ];
tol_pexpRB54s5=[10^-4.6, 10^-5.6 10^-8.6 10^-9.6];

for i=1:4
    [~,expHR53s3_sol,~,CPU,~,~] =expHR53s3(F,A,Jg,tspan(1),tspan(2),Y0,tol_expHR53s3(i),tol_expHR53s3(i),1e-10,100);
    expHR53s3_err(i)=norm(exac_sol-expHR53s3_sol,'inf');
    CPU_expHR53s3(i)=CPU;

    [~,expHR53s3a_sol,~,CPU,~,~]=expHR53s3a(F,A,Jg,tspan(1),tspan(2),Y0,tol_expHR53s3a(i),tol_expHR53s3a(i),1e-10,100);
    expHR53s3a_err(i)=norm(exac_sol-expHR53s3a_sol,'inf');
    CPU_expHR53s3a(i)=CPU;
    % 
    [~,expHR54s4_sol,~,CPU,~,~]=expHR54s4(F,A,Jg,tspan(1),tspan(2),Y0,tol_expHR54s4(i),tol_expHR54s4(i),1e-10,100);
    expHR54s4_err(i)=norm(exac_sol-expHR54s4_sol,'inf');
    CPU_expHR54s4(i)=CPU;

    [~,expHR54s4a_sol,~,CPU,~,~]=expHR54s4a(F,A,Jg,tspan(1),tspan(2),Y0,tol_expHR54s4a(i),tol_expHR54s4a(i),1e-10,100);
    expHR54s4a_err(i)=norm(exac_sol-expHR54s4a_sol,'inf');
    CPU_expHR54s4a(i)=CPU;

    [~,expHR54s5_sol,~,CPU,~,~]=expHR54s5(F,A,Jg,tspan(1),tspan(2),Y0,tol_expHR54s5(i),tol_expHR54s5(i),1e-10,100);
    expHR54s5_err(i)=norm(exac_sol-expHR54s5_sol,'inf');
    CPU_expHR54s5(i)=CPU;

    [~,expRB53s3_sol,~,CPU,~,~]=expRB53s3(F,A,g,Jg,tspan(1),tspan(2),Y0,tol_expRB53s3(i),tol_expRB53s3(i),1e-10,100);
    expRB53s3_err(i)=norm(exac_sol-expRB53s3_sol,'inf');
    CPU_expRB53s3(i)=CPU;

    [~,pexpRB54s4_sol,~,CPU,~,~]=pexpRB54s4(F,A,g,Jg,tspan(1),tspan(2),Y0,tol_pexpRB54s4(i),tol_pexpRB54s4(i),1e-10,100);
    pexpRB54s4_err(i)=norm(exac_sol-pexpRB54s4_sol,'inf');
    CPU_pexpRB54s4(i)=CPU;


    [~,pexpRB54s5_sol,~,CPU,~,~]=pexpRB54s5(F,A,g,Jg,tspan(1),tspan(2),Y0,tol_pexpRB54s5(i),tol_pexpRB54s5(i),1e-10,100);
    pexpRB54s5_err(i)=norm(exac_sol-pexpRB54s5_sol,'inf');
    CPU_pexpRB54s5(i)=CPU;

end

figure(1)
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 6.8 5.2]);
set(gcf,'PaperPositionMode','auto');

semilogy(CPU_expHR53s3, expHR53s3_err, 'v-','LineWidth',1, ...
    'Color',[0.5210 0.0860 0.8190],'MarkerSize',10); hold on   % purple down-triangle

semilogy(CPU_expHR53s3a, expHR53s3a_err, 'v-','LineWidth',1, ...
    'Color',[0.2310 0.6660 0.1960],'MarkerSize',10);           % green down-triangle

semilogy(CPU_expHR54s4, expHR54s4_err, '*-','LineWidth',1, ...
    'Color',[0.8660 0.3290 0.0000],'MarkerSize',10);           % orange star

semilogy(CPU_expHR54s4a, expHR54s4a_err, '*-','LineWidth',1, ...
    'Color',[0.9290 0.6940 0.1250],'MarkerSize',10);           % yellow star

semilogy(CPU_expHR54s5, expHR54s5_err, 'o-','LineWidth',1, ...
    'Color',[0.0660 0.4430 0.7450],'MarkerSize',10);           % blue circle

semilogy(CPU_expRB53s3, expRB53s3_err, '<-','LineWidth',1, ...
    'Color',[0.1840 0.7450 0.9370],'MarkerSize',10);           % cyan left-triangle

semilogy(CPU_pexpRB54s4, pexpRB54s4_err, '^-','LineWidth',1, ...
    'Color',[0.8190 0.0150 0.5450],'MarkerSize',10);           % magenta up-triangle

semilogy(CPU_pexpRB54s5, pexpRB54s5_err, '>-','LineWidth',1, ...
    'Color',[0.0660 0.4430 0.7450],'MarkerSize',10);           % blue right-triangle

xlabel('CPU time','FontSize',10)
ylabel('Error','FontSize',10)
grid on
box on

ax = gca;
ax.Units = 'normalized';
ax.Position = [0.24 0.16 0.72 0.72];
ax.FontSize = 10;
set(ax,'TickLength',2*get(ax,'TickLength'))

ymin = min([expHR53s3_err(:); expHR53s3a_err(:); expHR54s4_err(:); expHR54s4a_err(:); ...
            expHR54s5_err(:); expRB53s3_err(:); pexpRB54s4_err(:); pexpRB54s5_err(:)]);
ymax = max([expHR53s3_err(:); expHR53s3a_err(:); expHR54s4_err(:); expHR54s4a_err(:); ...
            expHR54s5_err(:); expRB53s3_err(:); pexpRB54s4_err(:); pexpRB54s5_err(:)]);

axis([9 19 ymin ymax])
set(ax,'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3])
axis square

lgd = legend('expHR53s3','expHR53s3a','expHR54s4','expHR54s4a', ...
             'expHR54s5','exprb53','pexprb54','pexprb55', ...
             'Location','westoutside');
lgd.FontSize = 10;

baseName = 'expHR_Schnakenberg_CPU5_embedded';
print(gcf,fullfile(outFolder, baseName), '-depsc');
saveas(gcf,fullfile(outFolder, [baseName '.png']));
addpath('../integrators','../phipmsimuliom','./Data')

function g = rhs(Y,nx,ny,a,b,kappa)

    N = nx*ny;

    U = Y(1:N);
    V = Y(N+1:end);

    U2V = (U.^2).*V;

    gu = kappa*(a - U + U2V);
    gv = kappa*(b - U2V);

    g = [gu; gv];

end

function J = Jac_g(Y,nx,ny,kappa)

    N = nx*ny;

    U = Y(1:N);
    V = Y(N+1:end);

    UV = U.*V;
    U2 = U.^2;

    I11 = [1,0;0,0];
    I12 = [0,1;0,0];
    I21 = [0,0;1,0];
    I22 = [0,0;0,1];

    J = kron(I11,spdiags(kappa*(-1 + 2*UV),0,N,N)) ...
      + kron(I12,spdiags(kappa*(U2),0,N,N)) ...
      + kron(I21,spdiags(-2*kappa*(UV),0,N,N)) ...
      + kron(I22,spdiags(-kappa*(U2),0,N,N));

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

