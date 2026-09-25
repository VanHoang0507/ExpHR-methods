function [t,expH_sol] = expH5s4(F,A,J,t0,t_end,u0,N,tol,c2,c3)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
    c4=(c3*5-3)/(10*c3-5);
   for i=1:N
        JnF=(A+J(u0))*F(u0);
        JF=J(u0)*F(u0);
        Un2=u0+c2*dt*F(u0)+c2^2*dt^2*phipm_simul_iom(1,dt*c2*A,[zero,zero,JnF],tol,1,2);
        Hn2=J(Un2)*F(Un2)-JF;
        Un3=u0+c3*dt*F(u0)+dt^2*phipm_simul_iom(1,dt*c3*A,[zero,zero,c3^2*JnF,c3^3/c2*Hn2],tol,1,2);%+dt^2*phipm_simul_iom(1,dt*c(2)*A,[zero,Hn2,Hn2],tol,1,2);
        Hn3=J(Un3)*F(Un3)-JF;
        Un4=u0+c4*dt*F(u0)+dt^2*phipm_simul_iom(1,dt*c4*A,[zero,zero,c4^2*JnF,c4^3/c2*Hn2],tol,1,2);%+dt^2*phipm_simul_iom(1,dt*c(2)*A,[zero,Hn2,Hn2],tol,1,2);
        Hn4=J(Un4)*F(Un4)-JF;
        u0=u0+dt*F(u0)+dt^2*phipm_simul_iom(1,dt*A,[zero,zero,JnF,c4*Hn3/(c3*c4-c3^2)+c3*Hn4/(c4*c3-c4^2),-2*Hn3/(c4*c3-c3^2)-2*Hn4/(c4*c3-c4^2)],tol,1,2);
   end
expH_sol=u0;
end
