
clearvars; close all; clc
%%% HW#1 problem 3 
%%% seasonal heat diffusion into a glacier
%%% Anna Valentine

% example pdepe solution
Ti = 0; % K
Th = 0; % K
Tc = 0; % K


%Ti = T0 + delT.*cos(omg*t); 


Dval = 1e-1; % diffusivity [length^2/time]

% definitions
x = linspace(-10,0,256); % mesh with L=1
t = linspace(0,10,256); 
sol = pdepe(0,@(x,t,u,DuDx)pdedef(x,t,u,DuDx,Dval),@(x)pdeic(x,Ti),@(xl,ul,xr,ur,t)pdebc(xl,ul,xr,ur,t,Th,Tc),x,t);

time = 5;
time_index = round(interp1(t,1:length(t),time));
plot(x,sol(time_index,:,1),'k','linewidth',2)
set(gca,'ticklabelinterpreter','latex','fontsize',18)

% % % % % % % % % % % % % % % % % % % % % % %
function [c,f,s] = pdedef(x,t,u,DuDx,Dval)
c = 1;
f = Dval*DuDx;
s = 0;
end
% % % % % % % % % % % % % % % % % % % % % % %
function u0 = pdeic(x,c0)
u0 = c0;
end
% % % % % % % % % % % % % % % % % % % % % % %
function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t,cl,cr)
T0 = 0;% K
delT = 1; % K
pl = ul-(T0 + delT.*cos(pi*t));
ql = 0;
pr = 0;
qr = -1;
end
% % % % % % % % % % % % % % % % % % % % % % %