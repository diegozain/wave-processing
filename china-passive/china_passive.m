% ------------------------------------------------
%
% apply interferometry and MASW 
% to passive seismic data from china.
%
% ------------------------------------------------ 


clear all
close all
clc

addpath('shared');

% % --------------------------------------------
% % load data 
% % --------------------------------------------
% 
% d1 = 'china-data/H0101.sac';
% d2 = 'china-data/H0102.sac';
% d3 = 'china-data/H0103.sac';
% d4 = 'china-data/H0104.sac';
% d5 = 'china-data/H0105.sac';
% d6 = 'china-data/H0106.sac';
% d7 = 'china-data/H0107.sac';
% d8 = 'china-data/H0108.sac';
% d9 = 'china-data/H0109.sac'; 
% d10 = 'china-data/H0110.sac';
% d11 = 'china-data/H0111.sac';
% d12 = 'china-data/H0112.sac'; 
% 
% [d1,hd] = rdSac(d1);
% [d2,hd] = rdSac(d2);
% [d3,hd] = rdSac(d3);
% [d4,hd] = rdSac(d4);
% [d5,hd] = rdSac(d5);
% [d6,hd] = rdSac(d6);
% [d7,hd] = rdSac(d7);
% [d8,hd] = rdSac(d8);
% [d9,hd] = rdSac(d9);
% [d10,hd] = rdSac(d10);
% [d11,hd] = rdSac(d11);
% [d12,hd] = rdSac(d12);
% 
% nt = numel(d1);
% d = zeros(nt,12);
% 
% d(:,1) = d1 ;
% d(:,2) = d2 ;
% d(:,3) = d3 ;
% d(:,4) = d4 ;
% d(:,5) = d5 ;
% d(:,6) = d6 ;
% d(:,7) = d7 ;
% d(:,8) = d8 ;
% d(:,9) = d9 ;
% d(:,10) = d10 ;
% d(:,11) = d11 ;
% d(:,12) = d12 ;
% 
% save('d.mat','d');

% --------------------------------------------
% load data
% --------------------------------------------
fprintf('\nload data \n')

d = load('d.mat');
d = d.d;

fs = 500; % [Hz]
dt = 1/fs; % [s]
dr = 10; % [m]

% --------------------------------------------
% process data
% --------------------------------------------
fprintf('\nbasic processing \n')

[nt,nr] = size(d);

dt = 1/fs;
fny = fs/2; % fny = 1/dt/2
T = (nt-1)*dt;
t = 0:dt:T;
df = fs/nt;
f = (0:nt-1) * df;
x = [0:nr-1]*dr;

d = normc(d);

% ---- see -----

figure;
hold on
for i=1:nr
plot(t,d(:,i)+i*0.005)
end
hold off
axis tight
xlabel('$t$ [s]' )
ylabel('receiver \#' )
title('$d(t,r)$ (raw. normalized.)' )
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
plot(f,d_pow,'.-')
xlabel('$f$ [Hz]' )
ylabel('$|d(f,r)|^2$ [m]' )
title('$|d(f,r)|^2$ before bandpass' )
fancy_figure()

% ------------------------------------------
%   bandpass
% ------------------------------------------
fprintf('\nbandpass \n')

f_low = 5e-0; % 1e-1
f_high = 2e+1; % 1e+2
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
plot(f,d_pow,'.-')
xlabel('$f$ [Hz]' )
ylabel('$|d(f,r)|^2$ [m]' )
title('$|d(f,r)|^2$ after bandpass' )
fancy_figure()

figure;
hold on
for i=1:nr
plot(t,d(:,i)+i*0.005)
end
hold off
axis tight
xlabel('$t$ [s]' )
ylabel('receiver \# [m]' )
title('$d(t,r)$ after bandpass' )
fancy_figure()

% figure;
% [clmap,c_axis,~] = fancy_colormap(d);
% imagesc(x,t,d)
% colormap(clmap)
% caxis(c_axis)
% ylabel('$t$ [s]' )
% xlabel('$r$ [m]' )
% title('$d(t,r)$ after bandpass' )
% fancy_figure()

% -------------------------------------------------------
% interferate
% -------------------------------------------------------
fprintf('\ninterferate \n')

d = normc(d);

% window data with ns 'sources'
% 
% choose ns from factor(nt)
% nt = 2^7  3^2  5^5
%
% want,
% T = nt_ * dt = dt*nt / ns
% 1/T = ns / dt*nt
% ns = dt*nt / T
%
ns = 720; % 2^3 * 3 * 5; % 5^5; % 2^3 * 3 * 5, 5^5 , 2^4 * 5^4 * 3;
d_windowed = window_d(d,ns);

% choose receiver to be virtual source
% index form (i.e. a # between 1 and nr)
%
virt = 12;

% virtual gather
%
% one receiver is virtual source
% stack over many sources
%
% d_windowed is a cube matrix (nt_, nr, ns)
% 		where nt_ is time of the ns sources.
%
[v_g, c_g, t_corr] = virt_gather(d_windowed,virt,dt);
c_g_ = sum(c_g,2);

% figure;
% [clmap,c_axis,~] = fancy_colormap(v_g);
% imagesc(x,t_corr,v_g)
% colormap(clmap)
% caxis(c_axis)
% xlabel('$r$ [m]' )
% ylabel('$t$ [s]' )
% title('$d_{\star}(t,r)$' )
% fancy_figure()

figure;
hold on
for i=1:nr
plot(t_corr,v_g(:,i)+i*0.3)
end
hold off
axis tight
xlabel('$t$ [s]' )
ylabel('receiver \#')
% ax = gca;
% ax.YTickLabel = {'1','2','3','4','5','6','7','8','9','10','11','12'};
title('$d_{\star}(t,r)$' )
fancy_figure()

% figure;
% subplot(2,3,[1,2,4,5])
% [clmap,c_axis,~] = fancy_colormap(c_g);
% imagesc(1:ns,t_corr,c_g)
% caxis(c_axis)
% colormap(clmap)
% xlabel('source No.' )
% ylabel('$t$ [s]' )
% title('correlation gather' )
% fancy_figure()
% 
% subplot(2,3,[3,6])
% plot(c_g_,t_corr,'k-')
% ax = gca;
% ax.YTick = [];
% axis tight
% title('pilot trace' )
% fancy_figure()

% -------------------------------------------------------
% get causal/acausal part
% -------------------------------------------------------

[v_g_causal,v_g_acausal,t_causal] = causal_acausal(v_g,t_corr);

v_g_causal = v_g_causal + v_g_acausal;
v_g_causal = normc(v_g_causal);

nt_causal = numel(t_causal);

% figure;
% [clmap,c_axis] = fancy_colormap(v_g_causal);
% imagesc(x,t_causal,v_g_causal)
% caxis(c_axis)
% colormap(clmap)
% xlabel('$r$ [m]' )
% ylabel('$t$ [s]' )
% title('$d_{\star}(t,r)$ causal' )
% fancy_figure()

figure;
hold on
for i=1:nr
plot(t_causal,v_g_causal(:,i)+i*0.05)
end
hold off
axis tight
xlabel('$t$ [s]' )
ylabel('receiver \#' )
% ax = gca;
% ax.YTickLabel = {'1','2','3','4','5','6','7','8','9','10','11','12'};
title('$d_{\star}(t,r)$ causal. before bandpass' )
fancy_figure()

% % -------------------------------------------------------
% % fk causal part
% % -------------------------------------------------------
% fprintf('\nfk analysis \n')
% 
% [d_fk,f,k,df,dk] = fk_rt(v_g_causal,dt,dr);
% d_fk = abs(d_fk);
% 
% figure;
% [clmap,~,amp] = fancy_colormap(d_fk);
% imagesc(k,f,d_fk)
% colormap(clmap)
% caxis([0 amp])
% xlabel('$k$ [1/m]' )
% ylabel('$f$ [Hz]' )
% title('$f-k$ of virtual gather' )
% fancy_figure()

% ------------------------------------------
%   bandpass
% ------------------------------------------
fprintf('\nbandpass \n')

v_g_causal = butter_butter(v_g_causal,dt,f_low,f_high,nbutter);

% taper edges to zero
% before fourier
%
for i=1:nr
  v_g_causal(:,i) = v_g_causal(:,i) .* tukeywin(nt_causal,0.001);
end
% d(t,r) -> d(f,r)
%
[v_g_causal_,f,df] = fourier_rt(v_g_causal,dt);
% power
%
v_g_causal_pow = abs(v_g_causal_).^2 / nt_causal;

v_g_causal = normc(v_g_causal);

figure;
plot(f,v_g_causal_pow,'.-')
ylabel('receiver \#' )
xlabel('$f$ [Hz]' )
title('$|d_{\star}(f,r)|^2$. after bandpass' )
fancy_figure()

figure;
hold on
for i=1:nr
plot(t_causal,v_g_causal(:,i)+i*0.1)
end
hold off
axis tight
xlabel('$t$ [s]' )
ylabel('receiver \#' )
% ax = gca;
% ax.YTickLabel = {'1','2','3','4','5','6','7','8','9','10','11','12'};
title('$d_{\star}(t,r)$ causal. after bandpass' )
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
f_max = 100;
vel_min = 100;
vel_max = 1000; % [m/s]
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
nt_ = fix(10*nt/T);
v_g_causal = v_g_causal(1:nt_,:);

[v_g_causal_,f,df] = fourier_rt(v_g_causal,dt);

% range of frequencies
%
f_disp = df:df:f_high;

% dispersion image ( nsx x nf_disp ) matrix
%
[disper_vxf,disper_sxf] = masw(v_g_causal_,x,sx,f,vx,f_disp);

% ---- see -----

figure;
[clmap,~,amp] = fancy_colormap(disper_vxf);
imagesc(f_disp,vx,disper_vxf)
colormap(clmap)
caxis([0 amp])
axis square
xlabel('$f$ [Hz]' )
ylabel('$v_x$' )
title('dispersion image' )
fancy_figure()


