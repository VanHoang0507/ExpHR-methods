function [t,expHR_sol,cpu]=expHR5s3(F,A,J,t0,t_end,u0,N,tol,c2)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
    % c2=1/3;
    c3=(c2*5-3)/(10*c2-5);
    b20=(c3/6-1/12)/(c2*(c3-c2));
    b30=(c2/6-1/12)/(c3*(c2-c3));
    tic
   for i=1:N
        Jn=(A+J(u0));
        U=phipm_simul_iom([c2 1]*dt,Jn,[zero,F(u0)],tol,1,2);
        Un2=u0+U(:,1);
        Hn2=(J(Un2)-J(u0))*F(Un2);
        varphi32=phipm_simul_iom(dt*c2,Jn,[zero,zero,zero,(b20/b30)*Hn2/(dt*c2)],tol,1,2);
        Un3=u0+phipm_simul_iom(dt*c3,Jn,[zero,F(u0),zero,(1/c2)*Hn2/dt],tol,1,2)+varphi32;
        Hn3=(J(Un3)-J(u0))*F(Un3);
        u0=u0+U(:,2)+phipm_simul_iom(dt,Jn,[zero,zero,zero,(c3*Hn2/(c2*c3-c2^2)+c2*Hn3/(c2*c3-c3^2))/dt,(-2*Hn2/(c2*c3-c2^2)-2*Hn3/(c2*c3-c3^2))/dt^2],tol,1,2);
   end
   cpu = toc;
expHR_sol=u0;
end
