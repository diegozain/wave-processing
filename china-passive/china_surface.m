% ------------------------------------------------
%
% d_vg are interferometric shot gathers 
% on the 12 receivers produced by 
% china_passive.m
% 
% purpose of code is to remove surface waves.
% not as easy as advertised. not working.
%
% ------------------------------------------------ 


clear all
close all
clc

addpath('shared');

% --------------------------------------------
% load data
% --------------------------------------------
fprintf('\nload data \n')

d = load('d_vg.mat');
% d = load('dx.mat');
d = d.d;

fs = 250; % [Hz]
dt = 1/fs; % [s]
dr = 10; % [m]

[nt,nr,ns] = size(d);

dt = 1/fs;
fny = fs/2; % fny = 1/dt/2
T = (nt-1)*dt;
t = 0:dt:T;
df = fs/nt;
f = (0:nt-1) * df;
r = [0:nr-1]*dr;

% ---- see -----

figure;
hold on
for i=1:nr
plot(t,d(:,i,5)+i*0.1)
end
hold off
axis tight
xlabel('$t$ [s]' )
ylabel('receiver \#' )
title('$d(t,r) \;\; \mbox{source} = 5$' )
fancy_figure()

figure;
[clmap,c_axis,~] = fancy_colormap(d(:,5:5+3,5));
imagesc( r(5:5+3), t, d(:,5:5+3,5));
caxis(c_axis);
colormap(clmap) 
xlabel('$r$ [m]'); 
ylabel('$t$ [s]');
title('$d(t,r) \;\; \mbox{source} = 5$');
fancy_figure()

% % -------------------------------------------------------
% % interferate
% % -------------------------------------------------------
% fprintf('\ninterferate \n')
% 
% % choose receiver to be virtual source
% % index form (i.e. a # between 1 and nr)
% %
% virt = 5;
% 
% % virtual gather
% %
% % one receiver is virtual source
% % stack over many sources
% %
% % d_windowed is a cube matrix (nt_, nr, ns)
% % 		where nt_ is time of the ns sources.
% %
% [v_g, c_g, t_corr] = virt_gather(d(:,:,1),virt,dt);
% c_g_ = sum(c_g,2);
% 
% figure;
% hold on
% for i=1:nr
% plot(t_corr,v_g(:,i)+i*0.3)
% end
% hold off
% axis tight
% xlabel('$t$ [s]' )
% ylabel('receiver \#')
% title('$d_{\star}(t,r). \;\; \mbox{source} = 1. \;\; r_{\bullet} = 5$' )
% fancy_figure()
% 
% figure;
% subplot(2,3,[1,2,4,5])
% [clmap,c_axis,~] = fancy_colormap(c_g);
% imagesc(1:ns,t_corr,c_g)
% caxis(c_axis)
% colormap(clmap)
% xlabel('source No.' )
% ylabel('$t$ [s]' )
% title('$\mbox{source} = 1. \;\; r_{\bullet} = 5$' )
% fancy_figure()
% 
% subplot(2,3,[3,6])
% plot(c_g_,t_corr,'k-')
% ax = gca;
% ax.YTick = [];
% axis tight
% title('pilot trace' )
% fancy_figure()
% 
% % -------------------------------------------------------
% % get causal/acausal part
% % -------------------------------------------------------
% 
% [v_g_causal,v_g_acausal,t_causal] = causal_acausal(v_g,t_corr);
% 
% v_g_causal = v_g_causal + v_g_acausal;
% v_g_causal = normc(v_g_causal);
% 
% nt_causal = numel(t_causal);
% 
% figure;
% hold on
% for i=1:nr
% plot(t_causal,v_g_causal(:,i)+i*0.1)
% end
% hold off
% axis tight
% xlabel('$t$ [s]' )
% ylabel('receiver \#' )
% title('$d_{\star}(t,r)$ causal. before bandpass' )
% fancy_figure()
% 
% % % ------------------------------------------
% % %   bandpass
% % % ------------------------------------------
% % fprintf('\nbandpass \n')
% % 
% % % v_g_causal = butter_butter(v_g_causal,dt,f_low,f_high,nbutter);
% % 
% % % taper edges to zero
% % % before fourier
% % %
% % for i=1:nr
% %   v_g_causal(:,i) = v_g_causal(:,i) .* tukeywin(nt_causal,0.001);
% % end
% % % d(t,r) -> d(f,r)
% % %
% % [v_g_causal_,f,df] = fourier_rt(v_g_causal,dt);
% % % power
% % %
% % v_g_causal_pow = abs(v_g_causal_).^2 / nt_causal;
% % 
% % v_g_causal = normc(v_g_causal);
% % 
% % figure;
% % plot(f,v_g_causal_pow,'.-')
% % ylabel('receiver \#' )
% % xlabel('$f$ [Hz]' )
% % title('$|d_{\star}(f,r)|^2$. after bandpass' )
% % fancy_figure()
% % 
% % figure;
% % hold on
% % for i=1:nr
% % plot(t_causal,v_g_causal(:,i)+i*0.1)
% % end
% % hold off
% % axis tight
% % xlabel('$t$ [s]' )
% % ylabel('receiver \#' )
% % title('$d_{\star}(t,r)$ causal. after bandpass' )
% % fancy_figure()
% 
% % -------------------------------------------------------
% % correct
% % -------------------------------------------------------
% 
% d_corrected = d(:,:,5) - v_g_causal;
% d_corrected = normc(d_corrected);
% 
% figure;
% hold on
% for i=1:nr
% plot(t_causal,d_corrected(:,i)+i*0.1)
% end
% hold off
% axis tight
% xlabel('$t$ [s]' )
% ylabel('receiver \#' )
% title('$d(t,r). \;\; \mbox{source} = 5.$ corrected' )
% fancy_figure()

% -------------------------------------------------------
%
% correct - all
%
% -------------------------------------------------------
fprintf('\ninterferate all\n')

% choose receiver to be virtual source
% index form (i.e. a # between 1 and nr)
%
virt = 5;
d_corrected = d(:,:,virt);

for i=[1:4 9:12]
  % virtual gather
  %
  % one receiver is virtual source
  % stack over many sources
  %
  % d_windowed is a cube matrix (nt_, nr, ns)
  % 		where nt_ is time of the ns sources.
  %
  [v_g, c_g, t_corr] = virt_gather(d(:,:,i),virt,dt);

  % -------------------------------------------------------
  % get causal/acausal part
  % -------------------------------------------------------
  [v_g_causal,v_g_acausal,t_causal] = causal_acausal(v_g,t_corr);
  v_g_causal = v_g_causal + v_g_acausal;
  v_g_causal = normc(v_g_causal);

  % -------------------------------------------------------
  % correct
  % -------------------------------------------------------
  d_corrected = d_corrected - v_g_causal;
  d_corrected = normc(d_corrected);
  d(:,:,virt) = d_corrected;
end

figure;
hold on
for i=1:nr
plot(t_causal,d_corrected(:,i)+i*0.1)
end
hold off
axis tight
xlabel('$t$ [s]' )
ylabel('receiver \#' )
title('$d(t,r). \;\; \mbox{source} = 5.$ corrected' )
fancy_figure()

figure;
[clmap,c_axis,~] = fancy_colormap(d_corrected(:,virt:virt+3));
imagesc( r(virt:virt+3), t, d_corrected(:,virt:virt+3));
caxis(c_axis);
colormap(clmap) 
xlabel('$r$ [m]'); 
ylabel('$t$ [s]');
title('$d(t,r) \;\; \mbox{source} = 5$. corrected');
fancy_figure()




