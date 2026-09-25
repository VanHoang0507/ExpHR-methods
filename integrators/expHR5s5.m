function [t,expHR5s5_sol,expHR5s5_cpu] = expHR5s5(F,A,J,t0,t_end,u0,N,tol,c2)    
dt = (t_end-t0)/N;
    t = linspace(t0, t_end, N+1);
    m = length(u0);
    zero=zeros(m,1);
    c3=0.3;
    % c3=1/4;
    % c4=2/3;
    % c5=1;
    c4=0.8;
    c5=1;
    d33=c4*c5/(c3*(c3-c4)*(c3-c5));
    d34=-2*(c4+c5)/(c3*(c3-c4)*(c3-c5));
    d35=6/(c3*(c3-c4)*(c3-c5));
    d43=c3*c5/(c4*(c4-c3)*(c4-c5));
    d44=-2*(c3+c5)/(c4*(c4-c3)*(c4-c5));
    d45=6/(c4*(c4-c3)*(c4-c5));
    d53=c3*c4/(c5*(c5-c3)*(c5-c4));
    d54=-2*(c3+c4)/(c5*(c5-c3)*(c5-c4));
    d55=6/(c5*(c5-c3)*(c5-c4));
    tic
   for i=1:N
        Ju0 = J(u0);
        Fu0 = F(u0);
        Jn=A+Ju0;
        Un2=u0+phipm_simul_iom(c2*dt,Jn,[zero,Fu0],tol,1,2);
        Hn2=(J(Un2)-Ju0)*F(Un2);
        U1=phipm_simul_iom([c3,c4,c5]*dt,Jn,[zero,F(u0),zero,1/c2*Hn2/dt],tol,1,2);
        Un3=u0+U1(:,1);
        Hn3=(J(Un3)-J(u0))*F(Un3);
        Un4=u0+U1(:,2);
        Hn4=(J(Un4)-J(u0))*F(Un4);
        Un5=u0+U1(:,3);
        Hn5=(J(Un5)-J(u0))*F(Un5);
        % c = [c3,c4,c5];
        % Un = zeros(m,3);
        % Hn = zeros(m,3);
        % maxNumCompThreads(1);
        % parfor j=1:3
        %     Un(:,j) = u0+phipm_simul_iom(c(j)*dt,Jn,[zero,Fu0,zero,1/c2*Hn2/dt],tol,1,2);
        %     F1 = F(Un(:,j));
        %     Hn(:,j)=(J(Un(:,j))*F1-Ju0*F1);
        % end
        % Hn3 = Hn(:,1);
        % Hn4 = Hn(:,2);
        % Hn5 = Hn(:,3);
        u0=u0+dt^2*phipm_simul_iom(1,dt*Jn,[zero,Fu0/dt,zero,d33*Hn3+d43*Hn4+d53*Hn5,d34*Hn3+d44*Hn4+d54*Hn5,d35*Hn3+d45*Hn4+d55*Hn5],tol,1,2);
   end
expHR5s5_cpu = toc;
expHR5s5_sol = u0;
end