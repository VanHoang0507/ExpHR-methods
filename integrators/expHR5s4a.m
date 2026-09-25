function [t,expHR_sol,cpu]=expHR5s4a(F,A,J,t0,t_end,u0,N,tol,c2)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
    c3 = c2;
    c4=(c3*5-3)/(10*c3-5);
    tic
   for i=1:N
        Ju0 = J(u0); 
        Fu0 = F(u0);
        Jn=A+Ju0;
        U=phipm_simul_iom([c2 c4]*dt,Jn,[zero,Fu0],tol,1,2);
        Un2 = u0 + U(:,1);
        Hn2=(J(Un2)-J(u0))*F(Un2);
        %U1 = phipm_simul_iom([c3,c4]*dt,Jn,[zero,F(u0),zero,1/c2*Hn2/dt],tol,1,2);
        Un3=u0+U(:,1)+dt^2*(c3^3/c2)/6*Hn2;
        Hn3=(J(Un3)-J(u0))*F(Un3);
        Un4=u0+U(:,2)+dt^2*(c4^3/c2)/6*Hn2;
        Hn4=(J(Un4)-J(u0))*F(Un4);
        u0=u0+dt^2*phipm_simul_iom(1,dt*Jn,[zero,Fu0/dt,zero,c4*Hn3/(c3*c4-c3^2)+c3*Hn4/(c4*c3-c4^2),-2*Hn3/(c4*c3-c3^2)-2*Hn4/(c4*c3-c4^2)],tol,1,2);
        % Jn=A+J(u0);
        % Un2=u0+phipm_simul_iom(0.4*dt,Jn,[zero,F(u0)],tol,1,2);
        % Hn2=(J(Un2)-J(u0))*F(Un2);
        % U1 = phipm_simul_iom([0.4,1]*dt,Jn,[zero,F(u0),zero,(1/0.4)*Hn2/dt],tol,1,2);
        % Un3=u0+U1(:,1);
        % Hn3=(J(Un3)-J(u0))*F(Un3);
        % Un4=u0+U1(:,2);
        % Hn4=(J(Un4)-J(u0))*F(Un4);
        % u0=u0+dt^2*phipm_simul_iom(1,dt*Jn,[zero,F(u0)/dt,zero,(25/6)*Hn3-(2/3)*Hn4,-(25/3)*Hn3+(10/3)*Hn4],tol,1,2);
   end
   cpu = toc;
expHR_sol=u0;
end
