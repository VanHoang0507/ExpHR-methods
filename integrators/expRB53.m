function [t,expRB53_sol,cpu]=expRB53(F,A,g,J,t0,t_end,u0,N,tol)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
tic
for i=1:N
      Jn=A+J(u0);
      gn=@(u) g(u)-J(u0)*u;
      U=phipm_simul_iom([1/2,9/10]*dt,Jn,[zero,F(u0)],tol,1,2);
      Un2=u0+U(:,1);
      Dn2=F(Un2)-Jn*Un2-gn(u0);
      varphi32=dt*phipm_simul_iom(1,0.5*dt*Jn,[zero,zero,zero,(27/25)*Dn2],tol,1,2);
      Un3=u0+U(:,2)+dt*phipm_simul_iom(1,(9/10)*dt*Jn,[zero,zero,zero,(729/125)*Dn2],tol,1,2)+varphi32;
      Dn3=F(Un3)-Jn*Un3-gn(u0);
      u0=u0+dt*phipm_simul_iom(1,dt*Jn,[zero,F(u0),zero,18*Dn2-(250/81)*Dn3,-60*Dn2+(500/27)*Dn3],tol,1,2);
end 
cpu = toc;
expRB53_sol=u0;
end