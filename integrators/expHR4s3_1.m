function [t,expHR4s3_sol,expHR4s3_cpu]=expHR4s3_1(F,A,J,t0,t_end,u0,N,tol,c)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
   tic 
   for i=1:N
        Jn=(A+J(u0));
        U=phipm_simul_iom([c(1),c(2)]*dt,Jn,[zero,F(u0)],tol,1,2);
        %U=phipm_simul_iom([1/4,1]*dt,Jn,[zero,F(u0)],tol,1,2);
        Un2=u0+U(:,1);
        Fun2 = F(Un2);
        Hn2=(J(Un2)-J(u0))*Fun2;
        Un3=u0+U(:,2);
        Hn3=(J(Un3)-J(u0))*F(Un3);
        u0=u0+dt^2*phipm_simul_iom(1,dt*Jn,[zero,F(u0),zero,(-c(1)*Hn3/(c(2)^2-c(1)*c(2))-c(2)*Hn2/(c(1)^2-c(1)*c(2))),(2*Hn3/(c(2)^2-c(1)*c(2))+2*Hn2/(c(1)^2-c(1)*c(2)))],tol,1,2);
   end
expHR4s3_cpu = toc;
expHR4s3_sol=u0;
end
