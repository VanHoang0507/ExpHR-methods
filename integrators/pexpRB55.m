function [t,expRB55_sol,expRB55_cpu]=pexpRB55(F,A,g,J,t0,t_end,u0,N,tol)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
tic
for i=1:N
      Jn=A+J(u0);
      gn=@(u) g(u)-J(u0)*u;
      U=phipm_simul_iom(1/2*dt,Jn,[zero,F(u0)],tol,1,2);
      Un2=u0+U;
      Dn2=F(Un2)-Jn*Un2-gn(u0);
      U1 = phipm_simul_iom([1/3,1/2,1]*dt,Jn,[zero,F(u0),zero,8*Dn2/dt^2],tol,1,2);
      Un3=u0+U1(:,2);
      Dn3=F(Un3)-Jn*Un3-gn(u0);
      Un4=u0+U1(:,1);
      Dn4=F(Un4)-Jn*Un4-gn(u0);
      Un5=u0+U1(:,3);
      Dn5=F(Un5)-Jn*Un5-gn(u0);
      % c = [1/3,1/2,1];
      % Un = zeros(m,3);
      % Dn = zeros(m,3);
      % maxNumCompThreads(1);
      % parfor j=1:3
      %     Un(:,j) = u0+phipm_simul_iom(c(j)*dt,Jn,[zero,F(u0),zero,8*Dn2/dt^2],tol,1,2);
      %     Dn(:,j)=F(Un(:,j))-Jn*Un(:,j)-gn(u0);
      % end
      % Dn3 = Dn(:,2);
      % Dn4 = Dn(:,1);
      % Dn5 = Dn(:,3);
      u0=u0+dt*phipm_simul_iom(1,dt*Jn,[zero,F(u0),zero,-32*Dn3+81*Dn4+Dn5,384*Dn3-729*Dn4-15*Dn5,-1152*Dn3+1944*Dn4+72*Dn5],tol,1,2);
end 
expRB55_cpu = toc;
expRB55_sol=u0;
end