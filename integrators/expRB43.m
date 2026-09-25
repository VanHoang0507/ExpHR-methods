function [t,expRB43_sol,CPU]=expRB43(F,A,g,J,t0,t_end,u0,N,tol)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
tic
for i=1:N
      Fu0 = F(u0);
      Ju0 = J(u0);
      Jn = (A+Ju0);
      gn =@(u) g(u)-Ju0*u;
      gnu0 = gn(u0);
      U=phipm_simul_iom([1/2,1]*dt,Jn,[zero,Fu0],tol,1,2);
      Un2=u0+U(:,1);
      Dn2=gn(Un2)-gnu0;
      Un3=u0+U(:,2)+phipm_simul_iom(dt,Jn,[zero,Dn2],tol,1,2);
      Dn3=gn(Un3)-gnu0;
      u0=u0+phipm_simul_iom(dt,Jn,[zero,Fu0,zero,(16*Dn2-2*Dn3)/dt^2,(-48*Dn2+12*Dn3)/dt^3],tol,1,2); 
end 
CPU = toc;
expRB43_sol=u0;
end

% Ju0 = J(u0);
%       Jnu0= A+Ju0;
%       gn=@(u) g(u)-J(u0)*u;
%       %Un2=u0+(2/4)*dt*phipm_simul_iom(1,(2/4)*dt*Jn(u0),[zero,F(u0)],tol,1,2);
%       %Un3=u0+dt*phipm_simul_iom(1,dt*Jn(u0),[zero,F(u0)],tol,1,2);
%       U=phipm_simul_iom([1/2,1]*dt,Jn,[zero,F(u0)],tol,1,2);
%       Un2=u0+U(:,1);
%       % Un3=u0+U(:,2);
%       Dn2=gn(Un2)-gn(u0);
%       Un3=u0+U(:,2)+phipm_simul_iom(dt,Jn,[zero,Dn2],tol,1,2);
%       Dn3=F(Un3)-Jn*Un3-gn(u0);
%       u0=u0+dt*phipm_simul_iom(1,dt*Jn,[zero,F(u0),zero,16*Dn2-2*Dn3,-48*Dn2+12*Dn3],tol,1,2); 