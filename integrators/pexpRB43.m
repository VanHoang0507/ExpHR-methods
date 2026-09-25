function [t,pexpRB43_sol,CPU]=pexpRB43(F,A,g,J,t0,t_end,u0,N,tol)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
tic
for i=1:N
      Fu0 = F(u0);
      Ju0 = J(u0);
      Jn = (A+Ju0);
      gn=@(u) g(u)-Ju0*u;
      gnu0 = gn(u0);
      U=phipm_simul_iom([1/2,1]*dt,Jn,[zero,Fu0],tol,1,2);
      Un2=u0+U(:,1);
      Dn2=F(Un2)-Jn*Un2-gnu0;
      Un3=u0+U(:,2);
      Dn3=F(Un3)-Jn*Un3-gn(u0);
      % c = [1/2,1];
      % Un = zeros(m,2);
      % Dn = zeros(m,2);
      % parfor j=1:2
      %     Un(:,j) = u0+phipm_simul_iom(c(j)*dt,Jn,[zero,F(u0)],tol,1,2);
      %     Dn(:,j) = F(Un(:,j))-Jn(Un(:,j))-gn(u0);
      % end
      % Dn2 = Dn(:,1);
      % Dn3 = Dn(:,2);
      u0=u0+phipm_simul_iom(dt,Jn,[zero,Fu0,zero,(16*Dn2-2*Dn3)/dt^2,(-48*Dn2+12*Dn3)/dt^3],tol,1,2); 
end 
CPU = toc;
pexpRB43_sol=u0;
end