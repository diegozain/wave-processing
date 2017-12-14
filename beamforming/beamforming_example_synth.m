% ------------------------------------------------
%
% apply beamforming to synthetic data
% on a linear (or modify d) array.
%
% ------------------------------------------------ 

close all
clear all

addpath('shared');

%-------------------
% synthetic data
%-------------------

% d is ( nt x nr )

% % wave coming from pi/2 all has arrivals at equal time.
% %
% d = zeros(length(r(:,1)),500);
% d(:,250) = ones(length(r(:,1)),1);

d = zeros(1100,8);
d(200,1) = 1;
d(300,2) = 1;
d(400,3) = 1;
d(500,4) = 1;
d(600,5) = 1;
d(700,6) = 1;
d(800,7) = 1;
d(900,8) = 1;

% % offset wave
% %
% d = eye(length(r(:,1)));
% d = [zeros(length(r(:,1)),2),d,zeros(length(r(:,1)),400)];

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
f_cutoff_high = 5;
[Ba,Bb] = butter(6,f_cutoff_high/fny,'low');
for i=1:nr
    d(:,i) = filtfilt(Ba,Bb,d(:,i));
end

% -------------------------------------------------------
% fourier: d(r,t), df -> d(r,f), f
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
r = [0 0; 1 0; 2 0; 3 0; 4 0; 5 0; 6 0; 7 0];

% see
%
figure;
hold on
plot(0,0,'.','MarkerSize',50)
plot(r(:,1),r(:,2),'.','MarkerSize',30)
hold off
grid on
axis tight
xlabel('x [Km]')
ylabel('y [Km]')
legend('origin','receivers')
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

% beamform this fucker
%
b = beamformer(fo,r,d_,f,sX,sY);
b = abs(b).^2 / nr;
% b = real(b);

figure;
imagesc(sx,sy,b)
axis image
xlabel('slowness x')
ylabel('slowness y')
title('beamformer')