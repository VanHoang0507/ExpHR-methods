function [t,expH_sol ] = expH4s31(F,A,J,t0,t_end,u0,N,tol,c2)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
    c3=(c2*5-5)/(10*c2-5);
    b20=(c3/6-1/12)/(c2*(c3-c2));
    b30=(c2/6-1/12)/(c3*(c2-c3));
   for i=1:N
        JnF=(A+J(u0))*F(u0);
        JF=J(u0)*F(u0);
        Un2=u0+c2*dt*F(u0)+c2^2*dt^2*phipm_simul_iom(1,dt*c2*A,[zero,zero,JnF],tol,1,2);
        Hn2=J(Un2)*F(Un2)-JF;
        varphi32=phipm_simul_iom(1,dt*c2*A,[zero,zero,zero,(c2^2*(b20/b30))*Hn2],tol,1,2);
        Un3=u0+c3*dt*F(u0)+dt^2*phipm_simul_iom(1,dt*c3*A,[zero,zero,c3^2*JnF,(c3^3/c2)*Hn2],tol,1,2)+dt^2*varphi32;
        Hn3=J(Un3)*F(Un3)-JF;
        u0=u0+dt*F(u0)+dt^2*phipm_simul_iom(1,dt*A,[zero,zero,JnF,c3*Hn2/(c2*c3-c2^2)+c2*Hn3/(c2*c3-c3^2),-2*Hn2/(c2*c3-c2^2)-2*Hn3/(c2*c3-c3^2)],tol,1,2);
   end
expH_sol=u0;
end
