function [t,expH_sol,cpu] = expH3(F,A,J,t0,t_end,u0,N,tol,c)
    dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
    tic
   for i=1:N
        JnF=(J(u0))*F(u0);
        Un2=u0+phipm_simul_iom(dt*c(1),A,[zero,F(u0),JnF],tol,1,2);
        Hn2=J(Un2)*F(Un2)-J(u0)*F(u0);
        u0=u0+phipm_simul_iom(dt,A,[zero,F(u0),JnF,(1/c(1))*Hn2/dt],tol,1,2); 
   end
   cpu = toc;
expH_sol=u0;
end
