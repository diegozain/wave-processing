% ------------------------------------------------
%
% thin data: from synthetic shot gather take
%            only traces from would be receivers.
%
% fatten data: from previous simulation, 
%              fatten data as if receivers had been 
%              everywhere.
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

d = load('data/electro/synth/dispersion/dispersion-1.mat');
% d = load('data/electro/synth/dispersion/dispersion-2.mat');
d = d.d_w_o;

[nt,nr] = size(d);

% output from readmala
%
T_ = 1.09334144418719048986e+02; % [ns]
%
% dt has to be multiplied by 1e-9 to enter wave-solver !!!!
%
dt_ = T_ / nt; 
fo = 4e+8; % [Hz] (mala would give in Ghz or some shit like that)
dr_ = 1.28762384550091570079e-02; % [m]

% derived from given data
%
fny = 1/dt_/2;
t_ = 0:dt_:T_-dt_;
df = 1/dt_/nt;
f = (0:nt-1) * df;
r_ = (0:nr-1) * dr_;

% --------------------------------------------
% thin data
% --------------------------------------------
fprintf('\nthin data \n')

% target parameters
%
epsilon_w = [7 2];

% known interval of target space
%
x = [0 20];
z = [x(1) 3];
t = T_;

% get x,z,t,r that will go into synthetic-wave.
%
[x,z,dx,dz,t,dt,r,dr] = w_xztr(epsilon_w,fo,x,z,t);

% -------
%
% this next step only has to be done for this particular 
% data set.
%
d = d(:,1:end-3);
%
% -------

% this part simulates what would happen in the field,
% real length 'x' (in nature) would be thinned by 'r'
% 
d_field = w_thin_d(d,x,r);

figure;

subplot(1,2,1)
imagesc(x,t_,d)
[clmap,c_axis,~] = fancy_colormap(d);
caxis(c_axis)
colormap(clmap)
c = colorbar;
c.TickLength = 0;
xlabel('$x$ [m]')
ylabel('$t$ [ns]')
title('$d(x,t)$')
fancy_figure()

subplot(1,2,2)
imagesc(r,t,d_field)
caxis(c_axis)
colormap(clmap)
c = colorbar;
c.TickLength = 0;
xlabel('$r$ [m]')
ylabel('$t$ [ns]')
title('$d(r,t)$')
fancy_figure()

% --------------------------------------------
% fatten data
% --------------------------------------------
fprintf('\nfatten data \n')

d_fat = w_fatten_d(d_field,x,r,t_);

figure;

subplot(1,2,1)
imagesc(x,t_,d_fat)
% [clmap,c_axis,~] = fancy_colormap(d_fat);
caxis(c_axis)
colormap(clmap)
c = colorbar;
c.TickLength = 0;
xlabel('$x$ [m]')
ylabel('$t$ [ns]')
title('$d(r,t)$ fattened to $d(x,t)$')
fancy_figure()

d_diff = d-d_fat;

subplot(1,2,2)
imagesc(x,t_,d_diff)
% [clmap,c_axis,~] = fancy_colormap(d_diff);
caxis(0.01*c_axis)
colormap(clmap)
c = colorbar;
c.TickLength = 0;
xlabel('$x$ [m]')
ylabel('$t$ [ns]')
title('$d(x,t)$ minus fattened $d(r,t)$')
fancy_figure()















