function [t,expHR43_sol,CPU]=expHR4s3(F,A,J,t0,t_end,u0,N,tol,c)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
   tic 
   for i=1:N
        Ju0 = J(u0);
        Fu0 = F(u0);
        Jn = (A+Ju0);
        U=phipm_simul_iom([c(1),c(2)]*dt,Jn,[zero,Fu0],tol,1,2);
        Un2=u0+U(:,1);
        Hn2=(J(Un2)*F(Un2)-Ju0*F(Un2));
        Un3=u0+U(:,2);
        Hn3=(J(Un3)*F(Un3)-Ju0*F(Un3));
        u0=u0+phipm_simul_iom(dt,Jn,[zero,Fu0,zero,(-c(1)*Hn3/(c(2)^2-c(1)*c(2))-c(2)*Hn2/(c(1)^2-c(1)*c(2)))/dt,(2*Hn3/(c(2)^2-c(1)*c(2))+2*Hn2/(c(1)^2-c(1)*c(2)))/dt^2],tol,1,2);
   end
   CPU = toc;
expHR43_sol=u0;
end
