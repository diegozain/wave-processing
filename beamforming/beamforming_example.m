% ------------------------------------------------
%
% apply beamforming to real data
% on a two-ring real array.
%
% ------------------------------------------------ 

close all
clear all

addpath('shared');

%----------------------------------------------------------
% real data
%----------------------------------------------------------

arr = load('data/array_data_example.mat');
d = arr.dta;
d = d';
% d is ( nt x nr ) matrix
[nt,nr] = size(d);

% sample rate is 200[Hz]
%
fs = 200;
fny = fs/2;
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

% -------------------------------------------------------
% fourier
% -------------------------------------------------------

[d_,f,df] = fourier_rt(d,dt);

% power
%
d_pow = abs(d_).^2 / nt;

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
title('data in time (detrend-demean-filter. cooked)')

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
% beamformer
% -------------------------------------------------------

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

% build slowness grid
%
vel_min = 1;
slow_max = 1/vel_min;
dsx = 0.01;
dsy = 0.01;

sx = -slow_max: dsx : slow_max;
sy = -slow_max: dsy : slow_max;

[sX,sY] = meshgrid(sx,sy);

% choose frequency
%
fo = 1;

% -----------------------
%
% beamform this fucker
%
% -----------------------

b = beamformer(fo,r,d_,f,sX,sY);
b = abs(b).^2 / nr;
% b = real(b);

% ---- see -----

figure;
imagesc(sx,sy,b)
axis image
xlabel('slowness x')
ylabel('slowness y')
title('beamformer')

% % plot like haney
% %
% wo = 2*pi*fo;
% beamPlot( b, sX, sY, slow_max, wo );


