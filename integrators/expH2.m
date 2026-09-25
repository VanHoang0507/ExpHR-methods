function [t,expH_sol]= expH2(F,A,J,t0,t_end,u0,N,tol)
    dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    n = length(u0);
    zero=zeros(n,1);
    for i=1:N
        JnF=J(u0)*F(u0);
        u0=u0+phipm_simul_iom(dt,A,[zero,F(u0),JnF],tol,1,2);
   end
expH_sol=u0;
end