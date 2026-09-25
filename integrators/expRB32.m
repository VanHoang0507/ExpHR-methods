function [t,expRB32_sol, CPU] = expRB32(F,A,g,J,t0,t_end,u0,N,tol)
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
      U = phipm_simul_iom(dt,Jn,[zero,Fu0],tol,1,2);
      Un2 = u0 + U(:,1);
      Dn2=F(Un2)-Jn*(Un2)-gnu0;
      u0=u0+U(:,1)+phipm_simul_iom(dt,Jn,[zero,zero,zero,2*Dn2/dt^2],tol,1,2); 
end  
CPU = toc;
expRB32_sol=u0;
end