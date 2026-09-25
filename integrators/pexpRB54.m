function [t,expRB54_sol,cpu]=pexpRB54(F,A,g,J,t0,t_end,u0,N,tol)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
tic
for i=1:N
      Jn=A+J(u0);
      gn=@(u) g(u)-J(u0)*u;
      Un2=u0+phipm_simul_iom((1/4)*dt,Jn,[zero,F(u0)],tol,1,2);
      Dn2=F(Un2)-Jn*Un2-gn(u0);
      U1 = phipm_simul_iom([1/2,9/10]*dt,Jn,[zero,F(u0),zero,32*Dn2/dt^2],tol,1,2);
      Un3=u0+U1(:,1);
      Dn3=F(Un3)-Jn*Un3-gn(u0);
      Un4=u0++U1(:,2);
      Dn4=F(Un4)-Jn*Un4-gn(u0);
      % c = [1/2,9/10];
      % Un = zeros(m,2);
      % Dn = zeros(m,2);
      % maxNumCompThreads(1);
      % parfor j=1:2
      %     Un(:,j) = u0+phipm_simul_iom(c(j)*dt,Jn,[zero,F(u0),zero,32*Dn2/dt^2],tol,1,2);
      %     Dn(:,j) = F(Un(:,j))-Jn*Un(:,j)-gn(u0);
      % end
      % Dn3 = Dn(:,1);
      % Dn4 = Dn(:,2);
      u0=u0+dt*phipm_simul_iom(1,dt*Jn,[zero,F(u0),zero,18*Dn3-(250/81)*Dn4,-60*Dn3+(500/27)*Dn4],tol,1,2);
end 
cpu =toc;
expRB54_sol=u0;
end