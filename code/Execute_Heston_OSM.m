clear all;clc;close all;clf;format long; format compact;

rho = 0.8; %% correlation factor
r = 0.03; %% interest rate
eta = 0.2; %% reversion level
kappa = 2; %% reversion rate
sig = 0.3; %% the volatility factor of volatility
T=1; %% maturity time
K=100; %% strike

ld = 0.5; %% Lambda

S=200; V=1;

% # of mesh grids
ns=40; nv=40; nt=40;
ds=S/ns; dv=V/nv; dt=T/nt;

% Test for blow-up in hybrid scheme
% nt = 15;

v=0:dv:V;
s=0:ds:S;
%%%%%%%%%%%%%%% Exact Solution at maturity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exc_u=max(0,s-K)'*ones(1,nv+1);

figure(1)
subplot(1,2,1)
mesh(0:ds:S,0:dv:V,exc_u')
axis tight
title('Exact solution at maturity')
xlabel('S')
ylabel('v')
zlabel('U')

%%%%%%%%%%%%%%% Exact Solution at initial time%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exc_iu = zeros(ns+1, nv+1);
for i = 1:ns+1
    percent = (i / (ns+1))*100;
    fprintf('Exact Loop : %4.1f %%\n', percent);
    for j = 1:nv+1
        exc_iu(i,j) = HestonCall(s(i),K,v(j),r,T,kappa,eta,sig,rho,ld);
    end
end

figure(1)
subplot(1,2,2)
mesh(0:ds:S,0:dv:V,exc_iu')
axis tight
title('Exact solution at initial time')
xlabel('S')
ylabel('v')
zlabel('U')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Heston model (Linear) %%%%%%%%%%%%%%%%%%%%%%%%
u = exc_u;
for t=1:nt
    percent = (t / nt)*100;
    fprintf('OSM Loop : %4.1f %%\n', percent);
    u_old=u; % Update u
    for j=2:nv
        sub_s = s(2:ns).^2*v(j)/2/ds^2-r*s(2:ns)/2/ds; % alpha
        main_s = -s(2:ns).^2*v(j)/ds^2-r/2-1/dt;       % beta
        sup_s = s(2:ns).^2*v(j)/2/ds^2+r*s(2:ns)/2/ds; % gamma
        
        c1=s(2).^2*v(j)/2/ds^2-r*s(2)/2/ds;             % alpha at s = 1
        cend=s(ns)^2*v(j)/2/ds^2+r*s(ns)/2/ds;          % gamma at s = Ns
        % LBC
        sup_s(1) = sup_s(1) - c1;
        main_s(1) = main_s(1) + 2*c1;
        sub_s(end) = sub_s(end) - cend;
        main_s(end) = main_s(end) + 2*cend;
        
        for i=2:ns
            b=-u_old((2:ns),j)'/dt; % bij
        end
        temp=Thomas(sub_s, main_s, sup_s, b);
        %u(i,j)=temp(i-1);
        u(2:ns,j) = temp(1:end);
        %end
        %%%%%%%%%%%%%%% Boundary %%%%%%%%%%%%%%%
        u(:,end)=2*u(:,end-1)-u(:,end-2);
        u(:,1)=2*u(:,2)-u(:,3);
        u(end,:)=2*u(end-1,:)-u(end-2,:);
        u(1,:)=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    u_old=u;
    for i=2:ns
        sub_v = sig^2*v(2:nv)/4/dv^2-kappa*(eta-v(2:nv))/2/dv; % alpha
        main_v = -sig^2*v(2:nv)/2/dv^2-r/2*ones(1,nv-1)-1/dt;  % beta
        sup_v = sig^2*v(2:nv)/4/dv^2+kappa*(eta-v(2:nv))/2/dv; % gamma
        
        c1t= sig^2*v(2)/4/dv^2-kappa*(eta-v(2))/2/dv;         % alpha at s = 1
        cendt = sig^2*v(nv)/4/dv^2+kappa*(eta-v(nv))/2/dv;    % gamma at s = Nv
        % LBC
        sup_v(1) = sup_v(1) - c1t;
        main_v(1) = main_v(1) + 2*c1t;
        sub_v(end) = sub_v(end) - cendt;
        main_v(end) = main_v(end) + 2*cendt;
        for j=2:nv
            b=-u_old(i,2:nv)/dt-rho*sig*s(i)*v(2:nv)/2*...
                (u_old(i+1,j+1)+u_old(i-1,j-1)-u_old(i+1,j-1)-u_old(i-1,j+1))...
                /(4*ds*dv);
        end
        temp=Thomas(sub_v, main_v, sup_v, b);
        u(i,2:nv) = temp(1:end);
        %   u(i,j)=temp(j-1);
        %end
        %%%%%%%%%%%%%%% Boundary %%%%%%%%%%%%%%%
        u(:,end)=2*u(:,end-1)-u(:,end-2);
        u(:,1)=2*u(:,2)-u(:,3);
        u(end,:)=2*u(end-1,:)-u(end-2,:);
        u(1,:)=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lin_u = u;
figure(2)
subplot(2,2,1)
mesh(0:ds:S,0:dv:V,lin_u')
axis tight
title('OSM method for Heston model with Linear BC')
xlabel('S')
ylabel('v')
zlabel('U')

figure(2)
subplot(2,2,3)
mesh(0:ds:S,0:dv:V,abs(lin_u - exc_iu)')
axis tight
title('Exact VS Linear BC')
xlabel('S')
ylabel('v')
zlabel('U')


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Heston model (New Hybrid) %%%%%%%%%%%%%%%%%%%%
clear u;
u = exc_u;
u_old = exc_u; % Update u
for t=1:nt
    percent = (t / nt)*100;
    fprintf('OSM Loop : %4.1f %%\n', percent);
    %%%%%%%%%%%%%%% Boundary %%%%%%%%%%%%%%%
    % Linear
    u_old(1, 2:end-1) = 2*u_old(2,2:end-1) - u_old(3,2:end-1);
    u_old(2:end-1, 1) = 2*u_old(2:end-1,2) - u_old(2:end-1,3);
    % Linear
    u_old(end,1:end-2) = 2*u_old(end-1,1:end-2) - u_old(end-2,1:end-2);
    u_old(1:end-2,end) = 2*u_old(1:end-2,end-1) - u_old(1:end-2,end-2);
    % New hybrid
    u_old(end,end) = 2*u_old(end-1,end-1) - u_old(end-2,end-2);
    u_old(end-1,end) = 2*u_old(end-2,end-1) - u_old(end-3,end-2);
    u_old(end,end-1) = 2*u_old(end-1,end-2) - u_old(end-2,end-3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=2:nv
        sub_s = s(2:ns).^2*v(j)/2/ds^2-r*s(2:ns)/2/ds; % alpha
        main_s = -s(2:ns).^2*v(j)/ds^2-r/2-1/dt;       % beta
        sup_s = s(2:ns).^2*v(j)/2/ds^2+r*s(2:ns)/2/ds; % gamma
        
        for i=2:ns
            b=-u_old((2:ns),j)'/dt; % bij
        end
        b(1) = b(1) - sub_s(1)*u_old(1,j);
        b(end) = b(end) - sup_s(end)*u_old(end,j);
        
        temp=Thomas(sub_s, main_s, sup_s, b);
        u(2:ns,j) = temp(1:end);
        %   u(i,j) = temp(i-1);
        %end
    end
    %%%%%%%%%%%%%%% Boundary %%%%%%%%%%%%%%%
    % Linear
    u(1, 2:end-1) = 2*u(2,2:end-1) - u(3,2:end-1);
    u(2:end-1, 1) = 2*u(2:end-1,2) - u(2:end-1,3);
    % Linear
    u(end,1:end-2) = 2*u(end-1,1:end-2) - u(end-2,1:end-2);
    u(1:end-2,end) = 2*u(1:end-2,end-1) - u(1:end-2,end-2);
    % New hybrid
    u(end,end) = 2*u(end-1,end-1) - u(end-2,end-2);
    u(end-1,end) = 2*u(end-2,end-1) - u(end-3,end-2);
    u(end,end-1) = 2*u(end-1,end-2) - u(end-2,end-3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=2:ns
        sub_v = sig^2*v(2:nv)/4/dv^2-kappa*(eta-v(2:nv))/2/dv; % alpha
        main_v = -sig^2*v(2:nv)/2/dv^2-r/2*ones(1,nv-1)-1/dt;  % beta
        sup_v = sig^2*v(2:nv)/4/dv^2+kappa*(eta-v(2:nv))/2/dv; % gamma
        
        for j=2:nv
            b=-u(i,2:nv)/dt-rho*sig*s(i)*v(2:nv)/2*...
                (u(i+1,j+1)+u(i-1,j-1)-u(i+1,j-1)-u(i-1,j+1))...
                /(4*ds*dv);
        end
        b(1) = b(1) - sub_v(1)*u(i,1);
        b(end) = b(end) - sup_v(end)*u(i,end);
        
        temp=Thomas(sub_v, main_v, sup_v, b);
        u_old(i,2:nv) = temp(1:end);
        %   u_old(i,j)=temp(j-1);
        %end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_u = u_old;
%%%%%%%%%%%%%%% Boundary %%%%%%%%%%%%%%%
% Linear
new_u(1, 2:end-1) = 2*new_u(2,2:end-1) - new_u(3,2:end-1);
new_u(2:end-1, 1) = 2*new_u(2:end-1,2) - new_u(2:end-1,3);
% Linear
new_u(end,1:end-2) = 2*new_u(end-1,1:end-2) - new_u(end-2,1:end-2);
new_u(1:end-2,end) = 2*new_u(1:end-2,end-1) - new_u(1:end-2,end-2);
% New hybrid
new_u(end,end) = 2*new_u(end-1,end-1) - new_u(end-2,end-2);
new_u(end-1,end) = 2*new_u(end-2,end-1) - new_u(end-3,end-2);
new_u(end,end-1) = 2*new_u(end-1,end-2) - new_u(end-2,end-3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(2)
subplot(2,2,2)
mesh(0:ds:S,0:dv:V,new_u')
axis tight
title('OSM method for Heston model with Hybrid BC')
xlabel('S')
ylabel('v')
zlabel('U')

figure(2)
subplot(2,2,4)
mesh(0:ds:S,0:dv:V,abs(new_u - exc_iu)')
axis tight
title('Exact VS New BC')
xlabel('S')
ylabel('v')
zlabel('U')

maxErrLin = max(max(abs(lin_u - exc_iu)))
rmseLin = sqrt(sum(sum(((lin_u - exc_iu).^2))) / ((nv+1)*(ns+1)))

maxErrNew = max(max(abs(new_u - exc_iu)))
rmseNew = sqrt(sum(sum(((new_u - exc_iu).^2))) / ((nv+1)*(ns+1)))