% ------------------------------------------------
%
% apply FTAN to synthetic shot gather.
%
% ------------------------------------------------ 

clear all
close all
clc

addpath('shared');

%-------------------
% synthetic data
%-------------------

% loads:
%
% D -> data
% dt
% offset -> receivers
load('data/fundamental_test_data.mat');

% d is ( nt x nr )
%
d = D;
[nt,nr] = size(d);
receivers = offset;
clear D offset

% sample rate is 200[Hz]
%
fs = 1/dt;
fny = fs/2;
T = (nt-1)*dt;
t = 0:dt:T;
df = fs/nt;
f = (0:nt-1) * df;

% ---- see -----

figure;
imagesc(receivers,t,d)
axis ij
[a,b]=caxis; 
caxis(0.1*[a,b]);
xlabel('receivers [m]')
ylabel('time [s]')
ylim([0 0.5]);
title('data in time (raw. like meat.)')

% -------------------------------------------------------
% fourier: d(r,t), df -> d(r,f), f
% -------------------------------------------------------

% [d_,f,df] = fourier_rt(d,dt);
% % [d_,f,df] = fourier_rt2(d,dt);
%
% % power
% %
% d_pow = abs(d_).^2 / nt;
%
% nfft = numel(f);
%
% f_ = makeFFTvector(dt,nfft);

% -------------------------------------------------------
% ftan
% -------------------------------------------------------

fmin = 1; % [Hz]
fmax = 100; % [Hz]
df = 2; % [Hz]

vmin = 100; % [m/s]
vmax = 350; % [m/s]
dv = 1; % [m/s]

% make velocity and frequency vectors
%
f = fmin : df : fmax;  % [Hz] frequencies to scan over
v = vmin : dv : vmax; % [m/s] velocities to consider

% choose width of gaussian ftan filter
%
alph = 5;

% choose receiver to process
%
dro = d(:,end);

% distance from source to receiver [km]
%
dsr = receivers(end);  

% ftan this guy
%
[disper_g_vf, disper_g_ft, dro_snr] = ftan(dro,dt,f,v,alph,dsr); 

% Plot ftan image
%
h = figure;
imagesc( f, v, disper_g_vf ); 
axis xy;
xlabel('frequency [Hz]'); 
ylabel('this should be group velocity [m/s]');
title('dispersion image. (f, gv)');

t = ( 0 : nt-1 ) .* dt;

% Plot ftan image
%
h = figure;
imagesc( f, t, abs(disper_g_ft)' );
axis xy;
ylabel('time [s]');
xlabel('frequency [Hz]');
title('dispersion image. (t, f)');

















