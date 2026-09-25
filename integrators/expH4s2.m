function [t,expH4s2_sol] = expH4s2(F,A,J,t0,t_end,u0,N,tol)
    dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
    c=0.5;
   for i=1:N
        JnF=(A+J(u0))*F(u0);
        Un2=u0+c(1)*dt*F(u0)+c(1)^2*dt^2*phipm_simul_iom(1,dt*c(1)*A,[zero,zero,JnF],tol,1,2);
        Hn2=J(Un2)*F(Un2)-J(u0)*F(u0);
        u0=u0+dt*F(u0)+dt^2*phipm_simul_iom(1,dt*A,[zero,zero,JnF,(1/c(1))*Hn2],tol,1,2);
   end
expH4s2_sol=u0;
end

