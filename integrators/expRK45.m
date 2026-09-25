function [t,expRK45_sol] = expRK45(F,A,g,t0,t_end,u0,N,tol)  
dt = (t_end-t0)/N;
t = linspace(t0, t_end, N+1);
m = length(u0);
zero=zeros(m,1);

for i=1:N
   Un2=u0+0.5*dt*phipm_simul_iom(1,dt*0.5*A,[zero,F(u0)],tol,1,2);
   Dn2=g(Un2)-g(u0);
   Un3=Un2+dt*phipm_simul_iom(1,dt*0.5*A,[zero,zero,Dn2],tol,1,2);
   Dn3=g(Un3)-g(u0);
   Un4=u0+dt*phipm_simul_iom(1,dt*A,[zero,F(u0),Dn2+Dn3],tol,1,2);
   Dn4=g(Un4)-g(u0);
   Un5=u0+dt*phipm_simul_iom(1,dt*0.5*A,[zero,0.5*F(u0),0.25*(2*Dn2+2*Dn3-Dn4),0.5*(-Dn2-Dn3+Dn4)],tol,1,2)+...
       dt*phipm_simul_iom(1,dt*A,[zero,zero,0.25*(Dn2+Dn3-Dn4),(-Dn2-Dn3+Dn4)],tol,1,2);
   Dn5=g(Un5)-g(u0);
   u0=u0+dt*phipm_simul_iom(1,dt*A,[zero,F(u0),-Dn4+4*Dn5,4*Dn4-8*Dn5],tol,1,2);
end 
  
expRK45_sol=u0;
end