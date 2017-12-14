% ------------------------------------------------
%
% apply MASW to synthetic data (electromagnetic).
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

dt = 0.02733; % [ ? ]
fs = 1/dt; % [ ? ]
dr = 0.012876; % [m]

% --------------------------------------------
% process data
% --------------------------------------------
fprintf('\nbasic processing \n')

[nt,nr] = size(d);

fny = fs/2; % fny = 1/dt/2
T = (nt-1)*dt;
t = 0:dt:T;
df = fs/nt;
f = (0:nt-1) * df;
x = [0:nr-1]*dr;

figure;
imagesc(x,t,normc(d))
colorbar
xlabel('$x$ [m]')
ylabel('$t$ [ns]')
title('data raw')
fancy_figure()

% d = basic_proce(d);

% ------------------------------------------
%   fourier
% ------------------------------------------
fprintf('\nfourier \n')

% taper edges to zero
% before fourier
%
for i=1:nr
  d(:,i) = d(:,i) .* tukeywin(nt,0.1);
end
% d(t,s) -> d(f,s)
%
[d_,f,df] = fourier_rt(d,dt);
% power
%
d_pow = abs(d_).^2 / nt;

figure;
hold on
for i=1:nr
plot(f,d_pow(:,i),'.-')
end
hold off
axis tight
xlabel('$f$')
ylabel('d power')
title('data in frequency. before bandpass')
fancy_figure()

% -------------------------------------------------------
% fk causal part
% -------------------------------------------------------
fprintf('\nfk analysis \n')

[d_fk,f,k,df,dk] = fk_rt(d,dt,dr);
d_fk = abs(d_fk);

figure;
imagesc(k,f,d_fk)
colormap(parula)
xlabel('$k$ [1/m]')
ylabel('$f$')
title('$f-k$ of data')
fancy_figure()

% -------------------------------------------------------
% MASW: d(x,fo) = sum( exp(...) * p(slow,fo) )_{slow_min,slow_max}
%
% is discretized as
%
% d = L p, L = L(x,slow,fo). -> find p given d_fo
% -------------------------------------------------------
fprintf('\nmasw \n')
% build array vector (has to be linear)
%
x_max = x(end);

% build slowness array
%
f_max = 0.3;
vel_min = 0.001;
vel_max = 0.1; % [m/s]
% dvx = 1e+1; % [m/s]

dsx   = 1 / x_max / f_max; % [s/m] compute minimum dsx
dsx   = dsx / 4; % go small enough not to alias

slow_max = 1 / vel_min; % [s/m]
slow_min = 1 / vel_max; % [s/m]

sx = slow_min : dsx : slow_max;
% vx = vel_min : dvx : vel_max;
vx = 1./sx;

% -------------------------------------------------------
% many freq
% -------------------------------------------------------

[d_,f,df] = fourier_rt(d,dt);

% power
%
d_pow = abs(d_).^2 / nt;

% figure;
% amp = max(v_g_pow(:));
% imagesc(x,f,v_g_pow,[-sc*0.01*amp sc*0.01*amp])
% colormap(clmap)
% ylabel('f [Hz]','Fontsize',20)
% xlabel('x [m]','Fontsize',20)
% title('virtual gather in frequency','Fontsize',20)

% range of frequencies
%
f_disp = df:df:0.3;

% dispersion image ( nsx x nf_disp ) matrix
%
[disper_vxf,disper_sxf] = masw(d_,x,sx,f,vx,f_disp);

% ---- see -----

figure;
[clmap,~,amp] = fancy_colormap(disper_vxf);
imagesc(f_disp,vx,disper_vxf)
colormap(clmap)
caxis([0 amp])
axis square
xlabel('$f$ [Hz]')
ylabel('$v_x$')
title('dispersion image')
fancy_figure()
