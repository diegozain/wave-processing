% ------------------------------------------------
%
% fatten data: from simulated field data that only
%            has traces where there were receivers,
%            fatten it as if there had been receivers
%            everywhere.
%
% ------------------------------------------------ 

close all
clear all
clc

addpath('shared');

% --------------------------------------------
% load data
% --------------------------------------------
fprintf('\nload data \n')

d = load('data/electro/synth/dispersion/dispersion-1-thin.mat');
d = d.d_field;

[nt,nr] = size(d);

% output from readmala
%
T = 1.09334144418719048986e+02; % [ns]
%
% dt has to be multiplied by 1e-9 to enter wave-solver !!!!
%
dt = T / nt; 
fo = 4e+8; % [Hz] (mala would give in Ghz or some shit like that)
dr = 6.43811922750457815701e-02; % [m]

% derived from given data
%
fny = 1/dt/2;
t = 0:dt:T-dt;
df = 1/dt/nt;
f = (0:nt-1) * df;
r = (0:nr-1) * dr;

% --------------------------------------------
% fatten data
% --------------------------------------------
fprintf('\nfatten data \n')

% target parameters
%
% real parameters for these data are 
% epsilon_w = [7 2]
%
epsilon_w = [18 2];

% known interval of target space
%
x = [0 1.9958169605e+01];
z = [x(1) 3];
t_ = T;

% get x,z that will go into synthetic-wave.
%
[x,z,dx,dz,t_,dt_,r_,dr_] = w_xztr(epsilon_w,fo,x,z,t_);

% fatten observed field data into wave-solver useful data
%
d_fat = w_fatten_d(d,x,r,t);

figure;

subplot(1,2,1)
imagesc(r,t,d)
[clmap,c_axis,~] = fancy_colormap(d);
caxis(c_axis)
colormap(clmap)
c = colorbar;
c.TickLength = 0;
xlabel('$r$ [m]')
ylabel('$t$ [ns]')
title('$d(r,t)$')
fancy_figure()

subplot(1,2,2)
imagesc(x,t,d_fat)
% [clmap,c_axis,~] = fancy_colormap(d_fat);
caxis(c_axis)
colormap(clmap)
c = colorbar;
c.TickLength = 0;
xlabel('$x$ [m]')
ylabel('$t$ [ns]')
title('$d(r,t)$ fattened to $d(x,t)$')
fancy_figure()















