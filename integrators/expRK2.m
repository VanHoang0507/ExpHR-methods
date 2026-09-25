function [t,expRK2_sol]=expRK2(F,A,g,t0,t_end,u0,N,tol,c2)
    dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    n = length(u0);
    zero=zeros(n,1);
    for i=1:N
        Un2=u0+c2*dt*phipm_simul_iom(1,c2*dt*A,F(u0),tol,1,2);
        Dn2=g(Un2)-g(u0);
        u0=u0+dt*phipm_simul_iom(1,dt*A,[zero,F(u0),(1/c2)*Dn2],tol,1,2);
    end
expRK2_sol=u0;
end
