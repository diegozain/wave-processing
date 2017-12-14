% ------------------------------------------------
%
% apply beamforming to real data
% on a two-ring real array,
% then perform interferometry.
%
% ------------------------------------------------ 


close all
clear all

addpath('shared');

%----------------------------------------------------------
% receivers
%----------------------------------------------------------

% receiver locations
%
receivers = load('data/stationCoordinatesYX_in_km.txt');
r = receivers;
r(:,1) = receivers(:,2);
r(:,2) = receivers(:,1);
nr = length(r(:,1));
clear receivers

% see
%
figure;
hold on
plot(0,0,'.k','MarkerSize',50)
for i=1:nr
plot(r(i,1),r(i,2),'.','MarkerSize',30)
end
hold off
grid on
axis([-1.2,1.2,-2,0.2])
xlabel('x [Km]')
ylabel('y [Km]')
legend({'origin','receivers'})
title('Array of receivers')
fancy_figure()

%----------------------------------------------------------
% real data
%----------------------------------------------------------

arr = load('data/array_data_example.mat');
d = arr.dta;
d = d';
[nt,nr] = size(d);
% d is ( nt x nr ) matrix

% sample rate is 200[Hz]
%
fs = 200;
fny = fs/2; % fny = 1/dt/2
dt = 1/fs;
T = (nt-1)*dt;
t = 0:dt:T;
df = fs/nt;
f = (0:nt-1) * df;

% ---- see -----

figure;
hold on
for i=1:nr
plot(t,d(:,i) + (i-1)*0.2)
end
hold off
axis tight
xlabel('time [s]')
ylabel('amplitude [a.u]')
title('data in time (raw. like meat.)')
fancy_figure()

%----------------------------------------------------------
% filter
%----------------------------------------------------------

d = basic_proce(d);

% high
%
% below 0.1 Hz is microseismic energy.
% but actually, nothing coherent comes
% out until 1 Hz.
%
f_cutoff_low = 1;
[Ba,Bb] = butter(6,f_cutoff_low/fny,'high');
for i=1:nr
    d(:,i) = filtfilt(Ba,Bb,d(:,i));
end

% low
%
% above 5 Hz is noise energy.
%
f_cutoff_high = 10;
[Ba,Bb] = butter(6,f_cutoff_high/fny,'low');
for i=1:nr
    d(:,i) = filtfilt(Ba,Bb,d(:,i));
end

% ---- see -----

figure;
hold on
for i=1:nr
plot(t,d(:,i) + (i-1)*0.2)
end
hold off
axis tight
xlabel('time [s]')
ylabel('amplitude [a.u]')
title('data in time (demean-detrend-filter. cooked)')
fancy_figure()

% -------------------------------------------------------
% interferate
% -------------------------------------------------------

% choose receiver to be virtual source
% index form (i.e. a # between 1 and nr)
%
virt = 1;

% virtual shot
%
% one receiver is virtual source
%
[v_r, t_corr] = virt_shot(d,virt,dt);

% ---- see -----

figure;
hold on
for i=1:nr
plot(t_corr,v_r(:,i) + (i-1)*10)
end
hold off
axis tight
xlabel('time [s]')
ylabel('amplitude [a.u]')
title('virtual shots')
fancy_figure()

% window data with 3 'sources'
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
ns = 3;
d_windowed = window_d(d,ns);
% T_ = 5;
% [d_windowed,ns] = window_dT(d,T_,dt);

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
% method 'interferate' uses xcorr2  | they ARE slightly
% method 'virt_gather' uses xcorr   | different (1e-16)
%
[v_g, c_g, t_corr] = interferate(d_windowed,virt,dt);
% [v_g, c_g, t_corr] = virt_gather(d_windowed,virt,dt);
c_g_ = sum(c_g,2);

figure;
[clmap,c_axis,~] = fancy_colormap(v_g);
imagesc(1:nr,t_corr,v_g)
colormap(clmap)
caxis(c_axis)
xlabel('receiver \#')
ylabel('time [s]')
title('virtual gather')
fancy_figure()

figure;

subplot(2,3,[1,2,4,5])
[clmap,c_axis,~] = fancy_colormap(c_g);
imagesc(1:ns,t_corr,c_g)
colormap(clmap)
caxis(c_axis)
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

% turn all receivers into virtual sources
%
nt_ = nt/ns;
ns_ = nr;
v_gs = zeros(2*nt_-1,nr,ns_);
for i=1:ns_
	virt = i;
  [v_g, ~, t_corr] = interferate(d_windowed,virt,dt);
	% [v_g, ~, t_corr] = virt_gather(d_windowed,virt,dt);
	v_gs(:,:,i) = v_g;
end

% virtual virtual gather
%
virt = nr-1;
[v_v_g, c_g_g, t_corr] = interferate(v_gs,virt,dt);
% [v_v_g, c_g_g, t_corr] = virt_gather(v_gs,virt,dt);

figure;
[clmap,c_axis,~] = fancy_colormap(v_v_g);
imagesc(1:nr,t_corr,v_v_g)
colormap(clmap)
caxis(c_axis)
% axis image
xlabel('receiver \#')
ylabel('time [s]')
title('virtual gather gather')
fancy_figure()

figure;
hold on
for i=nr-1:nr
plot(t_corr,v_v_g(:,i) + (i-1)*4e+4)
end
hold off
axis tight
xlabel('time [s]')
ylabel('amplitude [a.u]')
title('virtual gather gather. just two traces')
fancy_figure()
