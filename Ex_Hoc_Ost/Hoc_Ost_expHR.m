clc
close all
clear all
addpath('../integrators','../phipmsimuliom')
global m  X
%Space interval
a = 0; 

b = 1; 
%Time interval
t0 = 0; 

t_end = 1;

m = 400;

delta_x = (b-a)/(m);

x = linspace(a,b,m+1);

X=x(2:end-1)';

e = ones(m-1,1);
A = spdiags([e*1/delta_x^2 -2*e/delta_x^2 e*1/delta_x^2],[-1 0 1],m-1,m-1);
U0=X-X.^2;
Ju= @(U) spdiags(-2*U./(1+U.^2).^2,0,m-1,m-1);



A=blkdiag(0,A);



U0=[t0;U0]; %New U0

g = @(U) [1;(2+X-X.^2).*exp(U(1))-1./(1+(X.*(1-X).*exp(U(1))).^2)+1./(1+U(2:end).^2)];

F= @(U) A*U+g(U);

u_true = @(x,t) (x-x.^2)*exp(t) ;        

Jac=@(U)J(U);


N=[2,4,8,16,32];

for i=1:length(N)

    [t,expHR3_sol]=expHR3(F,A,Jac,0,1,U0,N(i),1e-10,0.45);
    expHR3_err(i)=norm(u_true(X,t(end))-expHR3_sol(2:end),'inf');

    [t,expHR4s2_sol]=expHR3(F,A,Jac,0,1,U0,N(i),1e-10,1/2);
    expHR4s2_err(i)=norm(u_true(X,t(end))-expHR4s2_sol(2:end),'inf');

    [t,expHR4s3_sol]=expHR4s3(F,A,Jac,0,1,U0,N(i),1e-10,[0.1,0.6]);
    expHR4s3_err(i)=norm(u_true(X,t(end))-expHR4s3_sol(2:end),'inf');

    [t,expHR5s3_sol]=expHR5s3(F,A,Jac,0,1,U0,N(i),1e-10,0.3);
    expHR5s3_err(i)=norm(u_true(X,t(end))-expHR5s3_sol(2:end),'inf');

    [t,expHR5s3a_sol]=expHR5s3a(F,A,Jac,0,1,U0,N(i),1e-10,0.3);
    expHR5s3a_err(i)=norm(u_true(X,t(end))-expHR5s3a_sol(2:end),'inf');

    [t,expHR5s4_sol]=expHR5s4(F,A,Jac,0,1,U0,N(i),1e-10,0.1);
    TDexpR5s4_err(i)=norm(u_true(X,t(end))-expHR5s4_sol(2:end),'inf');

    [t,expHR5s4a_sol]=expHR5s4a(F,A,Jac,0,1,U0,N(i),1e-10,0.4);
    expHR5s4a_err(i)=norm(u_true(X,t(end))-expHR5s4a_sol(2:end),'inf');

    [t,expHR5s5_sol]=expHR5s5(F,A,Jac,0,1,U0,N(i),1e-10,0.3); %(1/4,2/3,1 give order 6 )
    expHR5s5_err(i)=norm(u_true(X,t(end))-expHR5s5_sol(2:end),'inf');

end

expHR3_order=[];
expHR4s2_order=[];
expHR4s3_order=[];
expHR5s3_order=[];
expHR5s4_order=[];
expHR5s3a_order=[];
expHR5s4a_order=[];
expHR5s5_order=[];

for i=1:length(N)-1
    expHR3_order=[expHR3_order,log(expHR3_err(i)/expHR3_err(i+1))/log(2)];
    expHR4s2_order=[expHR4s2_order,log(expHR4s2_err(i)/expHR4s2_err(i+1))/log(2)];
    expHR4s3_order=[expHR4s3_order,log(expHR4s3_err(i)/expHR4s3_err(i+1))/log(2)];
    expHR5s3_order=[expHR5s3_order,log(expHR5s3_err(i)/expHR5s3_err(i+1))/log(2)];
    expHR5s4_order=[expHR5s4_order,log(TDexpR5s4_err(i)/TDexpR5s4_err(i+1))/log(2)];
    expHR5s3a_order=[expHR5s3a_order,log(expHR5s3a_err(i)/expHR5s3a_err(i+1))/log(2)];
    expHR5s4a_order=[expHR5s4a_order,log(expHR5s4a_err(i)/expHR5s4a_err(i+1))/log(2)];
    expHR5s5_order=[expHR5s5_order,log(expHR5s5_err(i)/expHR5s5_err(i+1))/log(2)];
end
expHR3_order
expHR4s2_order
expHR4s3_order
expHR5s3_order
expHR5s4_order
expHR5s3a_order
expHR5s4a_order
expHR5s5_order

   set(0,'DefaultTextFontSize',15)
   set(0,'DefaultAxesFontSize',15)
    h=1./N;
  h = 1./N;

figure(1)
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 6.8 5.2]);
set(gcf,'PaperPositionMode','auto');

loglog(h,expHR3_err,'v-','LineWidth',1,'MarkerSize',10); hold on
loglog(h,expHR4s2_err,'*-','LineWidth',1,'MarkerSize',10);
loglog(h,expHR4s3_err,'>-','LineWidth',1,'MarkerSize',10);
loglog(h,expHR5s3_err,'^-','LineWidth',1,'MarkerSize',10);
loglog(h,expHR5s3a_err,'^-','LineWidth',1,'MarkerSize',10);
loglog(h,TDexpR5s4_err,'d-','LineWidth',1,'MarkerSize',10);
loglog(h,expHR5s4a_err,'o-','LineWidth',1,'MarkerSize',9);
loglog(h,expHR5s5_err,'s-','LineWidth',1,'MarkerSize',10);

loglog(h,h.^3*0.1,'-+','LineWidth',1,'Color',"#77AC30");
loglog(h,h.^4*0.005,'-x','LineWidth',1,'Color',"#77AC30");
loglog(h,h.^5*0.0005,'--','LineWidth',1,'Color',"#77AC30");

xlabel('Number of time steps','FontSize',10)
ylabel('Error','FontSize',10)
grid on
box on

ax = gca;
ax.Units = 'normalized';
ax.Position = [0.24 0.16 0.72 0.72];
ax.FontSize = 10;
set(ax,'TickLength',2*get(ax,'TickLength'))
axis square
lgd = legend('expHR3s2','expHR4s2','expHR4s3', ...
        'expHR5s3','expHR5s3a','expHR5s4','expHR5s4a', ...
        'expHR5s5','Slope 3','Slope 4','Slope 5', ...
        'Location','westoutside');
    lgd.FontSize = 10;
set(ax,'XTick',fliplr(h))
xticklabels({'32','16','8','4','2'})

ymin = min([expHR3_err(:); expHR4s2_err(:); expHR4s3_err(:); ...
            expHR5s3_err(:); expHR5s3a_err(:); TDexpR5s4_err(:); ...
            expHR5s4a_err(:); expHR5s5_err(:)]);
ymax = max([expHR3_err(:); expHR4s2_err(:); expHR4s3_err(:); ...
            expHR5s3_err(:); expHR5s3a_err(:); TDexpR5s4_err(:); ...
            expHR5s4a_err(:); expHR5s5_err(:)]);

axis([1/32 1/2 ymin ymax])

set(ax,'YTick',[1e-13 1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2])

print('expHR_Order','-depsc');
rmpath('../integrators','../phipmsimuliom')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Jacobian]=J(U)
global m X
    Ju= @(U) spdiags(-2*U(2:end)./(1+U(2:end).^2).^2,0,m-1,m-1);
    K=Ju(U);
    Jacobian=blkdiag(0,K);
    Jacobian(2:end,1)=(2+X-X.^2).*exp(U(1))+2*(X.*(1-X).*exp(U(1))).^2./(1+(X.*(1-X).*exp(U(1))).^2).^2;

end

% function [t,TDexp_sol,cpu]=TDexpR5s3a(F,A,J,t0,t_end,u0,N,tol,c2)    
% dt = (t_end-t0)/N;
%     t = linspace(t0, t_end, N+1);
%     m = length(u0);
%     zero=zeros(m,1);
%     % c2=1/3;
%     c3=(c2*5-3)/(10*c2-5);
%     b20=(c3/6-1/12)/(c2*(c3-c2));
%     b30=(c2/6-1/12)/(c3*(c2-c3));
%     tic
%    for i=1:N
%         Jn=(A+J(u0));
%         U = phipm_simul_iom([c2, c3]*dt,Jn,[zero,F(u0)],tol,1,2);
% 
%         Un2 = u0 + U(:,1);
%         Hn2=(J(Un2)-J(u0))*F(Un2);
% 
%         Un3 = u0 + U(:,2)+ dt^2*(c2^2*b20/(6*b30)+c3^3/(6*c2))*Hn2;
%         Hn3=(J(Un3)-J(u0))*F(Un3);
%         u0=u0+phipm_simul_iom(dt,Jn,[zero,F(u0),zero,(c3*Hn2/(c2*c3-c2^2)+c2*Hn3/(c2*c3-c3^2))/dt,(-2*Hn2/(c2*c3-c2^2)-2*Hn3/(c2*c3-c3^2))/dt^2],tol,1,2);
%    end
%    cpu = toc;
% TDexp_sol=u0;
% end
