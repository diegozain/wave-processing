% ------------------------------------------------
%
% apply MASW to synthetic data (seismic).
%
% ------------------------------------------------ 

close all
clear all

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

[d_,f,df] = fourier_rt(d,dt);
% [d_,f,df] = fourier_rt2(d,dt);

% power
%
d_pow = abs(d_).^2 / nt;

% ---- see -----

figure;
hold on
for i=1:nr
plot(f,d_pow(:,i),'.-')
end
hold off
axis tight
xlabel('f')
ylabel('d power')
title('data in frequency')

% -------------------------------------------------------
% MASW: d(x,fo) = sum( exp(...) * p(slow,fo) )_{slow_min,slow_max}
%
% is discretized as
%
% d = L p, L = L(x,slow,fo). -> find p given d_fo
% -------------------------------------------------------

% build array vector (has to be linear)
%
x = receivers;
x_max = x(end);

% build slowness array
%
f_max = 100;
vel_min = 200;
vel_max = 1000; % [m/s]
dvx = 1; % [m/s]

dsx   = 1 / x_max / f_max; % [s/m] compute minimum dp
dsx   = dsx / 4; % go small enough not to aliase

slow_max = 1 / vel_min; % [s/m]
slow_min = 1 / vel_max; % [s/m]

sx = slow_min : dsx : slow_max;
% vx = vel_min : dvx : vel_max;
vx = 1./sx;

% -------------------------------------------------------
% many freq
% -------------------------------------------------------

% range of frequencies
%
f_disp = df:df:f_max+df;

% dispersion image ( nsx x nf_disp ) matrix
%
[disper_vxf,disper_sxf] = masw(d_,x,sx,f,vx,f_disp);

% ---- see -----

figure;
imagesc(f_disp,vx,disper_vxf)
colorbar
axis square
xlabel('frequency')
ylabel('veolocity x')
title('dispersion image')

figure;
imagesc(f_disp,sx,disper_sxf)
colorbar
axis square
xlabel('frequency')
ylabel('slowness x')
title('dispersion image')
