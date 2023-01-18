clearvars; close all; clc

% example pdepe solution
Ti = 100; % K
T0 = 0; % K
TL = 0; % K

K = 1e-1; % diffusivity [length^2/time]
L = 15;    % L distance [length]
x = linspace(0,L,256); % mesh with L=15
t = linspace(0,180,256); % time = 10 s 
s = .1;   % we have a source term (uniform heating) [degree/time]
%n = linspace(1,100,100); % steps of Cm 
tA = 180; % I think time for the analytical does not need to be a step 

%%%%%%% ANAYLYTICAL SOLN %%%%%

%%%%% for t = 10 s, T = ??
T = T_analytical(K, L, s, t, x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOLN
sol = pdepe(0,@(x,t,u,DuDx)pdedef(x,t,u,DuDx,K),@(x)pdeic(x,Ti),@(xl,ul,xr,ur,t)pdebc(xl,ul,xr,ur,t,T0,TL),x,t);

time = 180;
time_index = round(interp1(t,1:length(t),time));
plot(x,sol(time_index,:,1),'k','linewidth',2)
hold on
plot(x,(T*(.1)),'b','linewidth',2)   %%% My analytical solution doesn't work but it's kinda close ! 
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',18)
legend('.pdepe soln.','analytical')
title('Q1: too diffuse')
xlabel('x-position') 
ylabel('Temperature (K)') 

% % % % % % % % % % % % % % % % % % % % % % %
function [c,f,s] = pdedef(x,t,u,DuDx,K)
c = 1;
f = K*DuDx;
s = 1;   % we have a source term (uniform heating) [degree/time]
end
% % % % % % % % % % % % % % % % % % % % % % %
function u0 = pdeic(x,c0)
u0 = c0;
end
% % % % % % % % % % % % % % % % % % % % % % %
function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t,cl,cr)
pl = ul-cl;
ql = 0;
pr = ur-cr;
qr = 0;
end
% % % % % % % % % % % % % % % % % % % % % % %
% compute Cn, sum in loop for each Cn, find T for each x using t. 
function Cn = cM_fun(L, K, s, n)
    term1= 1./((pi.^3).*(n.^3)) ;
    term2 = L.*cos(pi.*n)-1;
    term3 = (s.*(L.^2)./K)-((pi.^2).*(n.^2)); 
    term4 = (pi.*s.*(L.^3).*n.*(1./K).*sin(n.*pi)); 
    Cn = (2/L)*(term1.*term2.*term3)+term4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute T(x,t) 
function T = T_analytical(K, L, s, t, x)
    T=0;
    for n=1:150
        Cn = cM_fun(L, K, s, n);
        term1 = exp(-(n.^2).*(pi.^2).*K.*t.*(1./(L.^2)));
        term2 = sin(n.*pi.*x.*(1./L));
        term3 = (-s.*(x.^2).*(1./(2.*K)));
        term4 = (s.*L.*x)./(2.*K);
        T = T + (Cn.*term1.*term2)+term3+term4; 
    end

end
    