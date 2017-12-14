% ------------------------------------------------
%
% apply interferometry and FTAN to 
% passive seismic data.
%
% ------------------------------------------------ 

clear all
close all
clc

addpath('shared')
addpath('ERB_data')

% % --------------------------------------------
% % load data & correct gaps and overlaps
% % --------------------------------------------
% 
% % Stations are
% % BSM1 (center),
% % BSM3 (west),
% % BSM7 (east, alumni center station)
% 
% % loads data, dt and fs.
% %
% load_data;
% 
% % --------------------------------------------
% % bin selected traces
% % --------------------------------------------
% 
% % gives data d to use before filtering
% %
% bin_data;
% save('d.mat','d');

% --------------------------------------------
% load data
% --------------------------------------------
fprintf('\nload data \n')

d = load('dz.mat');
% d = load('dx.mat');
d = d.d;

fs = 250; % [Hz]
dt = 1/fs; % [s]
dr = 177; % [m]

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
r = [0:nr-1]*dr;

figure;
hold on
plot(t,d(:,1),'b-')
plot(t,d(:,2),'r-')
hold off
axis tight
xlabel('$t$ [s]')
ylabel('counts')
legend({'west','east'},'Location','best')
title('$d(t,r)$ raw')
fancy_figure()

d = basic_proce(d);

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
xlabel('$f$ [Hz]')
ylabel('$d$ power')
legend({'west','east'},'Location','best')
title('$|d(f,r)|^2$ before bandpass')
fancy_figure()

% ------------------------------------------
%   bandpass
% ------------------------------------------
fprintf('\nbandpass \n')

f_low = 1e-0; % 1e-1
f_high = 1e+2; % 1e+2 5e+0
nbutter = 6;
d = butter_butter(d,dt,f_low,f_high,nbutter);

% taper edges to zero
% before fourier
%
for i=1:nr
  d(:,i) = d(:,i) .* tukeywin(nt,0.1);
end
% d(t,r) -> d(f,r)
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
xlabel('$f$ [Hz]')
ylabel('$d$ power')
legend({'west','east'},'Location','best')
title('$|d(f,r)|^2$ after bandpass')
fancy_figure()

figure;
hold on
plot(t,d(:,1),'b-')
plot(t,d(:,2),'r-')
hold off
axis tight
xlabel('$t$ [s]')
ylabel('amplitude [a.u]')
legend({'west','east'},'Location','best')
title('$d(t,r)$ after bandpass')
fancy_figure()

% -------------------------------------------------------
% interferate
% -------------------------------------------------------
fprintf('\ninterferate \n')

d = normc(d);

% window data with ns 'sources'
%
% choose ns from factor(nt)
%
% want,
% T_ = nt_ * dt = dt*nt / ns
% 1/T_ = ns / dt*nt
% ns = dt*nt / T_
%
% function 'window_d' uses time-index as windowing criteriea,
% not time.
%
% function 'window_dT' uses time as windowing criteriea,
% not time-index.
%
% ns = 257;
% d_windowed = window_d(d,ns);
T_ = 10;
[d_windowed,ns] = window_dT(d,T_,dt);

% choose receiver to be virtual source
% index form (i.e. a # between 1 and nr)
%
virt = 1;

% virtual gather
%
% one receiver is virtual source
% stack over many sources
%
% d_windowed is a cube matrix (nt_, nr, ns)
% 		where nt_ is time of the ns sources.
% 
% method 'interferate' uses xcorr2  - they ARE slightly
% method 'virt_gather' uses xcorr   - different (1e-16)
%
[v_g, c_g, t_corr] = interferate(d_windowed,virt,dt);
c_g_ = sum(c_g,2);

figure;

subplot(2,2,[1,2])
plot(t_corr,v_g(:,1),'b-')
axis tight
ylabel('amplitude [a.u]')
legend({'west'},'Location','best')
title('virtual gather $d_{\star}(t,r)$')
fancy_figure()

subplot(2,2,[3,4])
plot(t_corr,v_g(:,2),'r-')
axis tight
ylabel('amplitude [a.u]')
xlabel('$t$ [s]')
legend({'east'},'Location','best')
fancy_figure()



figure;

subplot(2,3,[1,2,4,5])
[clmap,c_axis,~] = fancy_colormap(c_g);
imagesc(1:ns,t_corr,c_g)
colormap(clmap)
caxis(0.01*c_axis)
c = colorbar;
c.TickLength = 0;
xlabel('source \#')
ylabel('$t$ [s]')
title('correlation gather')
fancy_figure()

subplot(2,3,[3,6])
plot(c_g_,t_corr,'k-')
axis tight
ax = gca;
ax.YTick = [];
title('pilot trace')
fancy_figure()

% -------------------------------------------------------
% get causal/acausal part
% -------------------------------------------------------

[v_g_causal,v_g_acausal,t_causal] = causal_acausal(v_g,t_corr);

v_g_causal = v_g_causal + v_g_acausal;
v_g_causal = normc(v_g_causal);

nt_causal = numel(t_causal);

figure;

subplot(2,2,[1,2])
plot(t_causal,v_g_causal(:,1),'b-')
axis tight
ylabel('amplitude [a.u]')
legend({'west'},'Location','best')
title('$d_{\star}(t,r)$ causal part (folded)')
fancy_figure()

subplot(2,2,[3,4])
plot(t_causal,v_g_causal(:,2),'r-')
axis tight
ylabel('amplitude [a.u]')
xlabel('$t$ [s]')
legend({'east'},'Location','best')
fancy_figure()

% ------------------------------------------
%   bandpass
% ------------------------------------------
fprintf('\nbandpass \n')

f_low = 1e-0; % 1e-1
f_high = 1e+1; % 1e+2 5e+0
nbutter = 6;
v_g_causal = butter_butter(v_g_causal,dt,f_low,f_high,nbutter);

% taper edges to zero
% before fourier
%
for i=1:nr
  v_g_causal(:,i) = v_g_causal(:,i) .* tukeywin(nt_causal,0.1);
end
% d(t,r) -> d(f,r)
%
[v_g_causal_,f,df] = fourier_rt(v_g_causal,dt);
% power
%
v_g_causal_pow = abs(v_g_causal_).^2 / nt_causal;

v_g_causal = normc(v_g_causal);
v_g_causal_pow = normc(v_g_causal_pow);

figure;
hold on
for i=1:nr
plot(f,v_g_causal_pow(:,i),'.-')
end
hold off
axis tight
xlabel('$f$ [Hz]')
ylabel('$d_{\star}$ power')
legend({'west','east'},'Location','best')
title('$|d_{\star}(f,r)|^2$ after bandpass')
fancy_figure()

figure;

subplot(2,2,[1,2])
plot(t_causal,v_g_causal(:,1),'b-')
axis tight
ylabel('amplitude [a.u]')
legend({'west'},'Location','best')
title('$d_{\star}(t,r)$ after bandpass')
fancy_figure()

subplot(2,2,[3,4])
plot(t_causal,v_g_causal(:,2),'r-')
axis tight
ylabel('amplitude [a.u]')
xlabel('$t$ [s]')
legend({'east'},'Location','best')
fancy_figure()

% -------------------------------------------------------
% ftan
% -------------------------------------------------------
fprintf('\nftan \n')

fmin = 2;  % [Hz] 25 2
fmax = 12; % [Hz] 80 12
% df = 0.5;  % [Hz] 1 0.5

vmin = 1e+2; % [m/s]
vmax = 1e+3; % [m/s]
dv = 1e+1;   % [m/s]

% make velocity and frequency vectors
%
f = fmin : df : fmax;  % [Hz] frequencies to scan over
v = vmin : dv : vmax;  % [m/s] velocities to consider

% choose width of gaussian ftan filter
%
alph = 5;

% choose receiver to process
%
dro = v_g_causal(:,2);

% distance from source to receiver [m]
%
dsr = 177;  

% ftan this guy
%
[disper_g_vf, disper_g_ft, dro_snr] = ftan(dro,dt,f,v,alph,dsr); 

% Plot ftan image
%
figure;
[clmap,~,amp] = fancy_colormap(disper_g_vf);
imagesc( f, v, disper_g_vf);
caxis([0 amp])
colormap(clmap) 
axis xy;
xlabel('$f$ [Hz]'); 
ylabel('$v_g$ [m/s]');
title('dispersion image');
fancy_figure()





