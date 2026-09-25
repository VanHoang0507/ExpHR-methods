function [t,expRB42_sol, CPU]=expRB42(F,A,g,J,t0,t_end,u0,N,tol)
    dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
tic
for i=1:N
      Fu0 = F(u0);
      Ju0 = J(u0); 
      Jn =  (A+Ju0);
      gn =@(u) g(u)-Ju0*u;
      gnu0 = gn(u0);
      U = phipm_simul_iom([(3/4) 1]*dt,Jn,[zero,Fu0],tol,1,2);
      Un2 = u0 + U(:,1);
      Dn2=F(Un2)-Jn*(Un2)-gnu0;
      u0=u0+U(:,2)+phipm_simul_iom(dt,Jn,[zero,zero,zero,32/9*Dn2/dt^2],tol,1,2); 
end  
CPU = toc;
expRB42_sol=u0;
end

% Jn=@(u) A+J(u);
%       gn=@(u) g(u)-J(u0)*u;
%       Un2=u0+(3/4)*dt*phipm_simul_iom(1,(3/4)*dt*Jn(u0),[zero,F(u0)],tol,1,2);
%       Dn2=F(Un2)-Jn(u0)*Un2-gn(u0);
%       u0=u0+dt*phipm_simul_iom(1,dt*Jn(u0),[zero,F(u0),zero,32/9*Dn2],tol,1,2); 