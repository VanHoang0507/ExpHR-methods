function [t,expHR_sol,CPU]=expHR3(F,A,J,t0,t_end,u0,N,tol,c)
    dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
   tic
   for i=1:N
        Ju0 = J(u0);
        Jn = (A+Ju0);
        Fu0 = F(u0);
        % U = phipm_simul_iom(dt,Jn,[zero,Fu0],tol,1,2);
        % Un2=u0+U;
        % Hn2=(J(Un2)*F(Un2)-Ju0*F(Un2));
        % u0=u0+U+phipm_simul_iom(dt,Jn,[zero,zero,zero,Hn2/dt],tol,1,2);

   
        Un2=u0+phipm_simul_iom(c*dt,Jn,[zero,Fu0],tol,1,2);
        Fun2 = F(Un2);
        Hn2=(J(Un2)-Ju0)*Fun2;
        u0=u0++phipm_simul_iom(dt,Jn,[zero,Fu0,zero,(1/c)*Hn2/dt],tol,1,2);
   end
   CPU = toc;
expHR_sol=u0;
end
