% ------------------------------------------------
%
% apply preprocessing and 
% image in time with gain-bandpass 
% a gpr radargram.
%
% ------------------------------------------------ 

close all
clear all

addpath('shared');

% for plotting
%
sc=0.3;
clo=[0 0 1];cmid=[.99 .99 .99];chi=[.5 .12 0];
clmap=mkclrmap(clo,cmid,chi);

% ------------------------------------------
%   read raw data
% ------------------------------------------

[d,dt,f0,dr] = readmala('data/electro/benin/Profile5');
% dt = 0.4; % it actually is dt = 0.38
% f0 = 200; % it actually is f0 = 250

[nt,nr] = size(d);
r = [0:nr-1]*dr;
t = [0:nt-1]*dt;
fny = 1/dt/2;

figure;
amp = max(max(d));
imagesc(r,t,d,[-sc*amp sc*amp])
colormap(clmap)
ylabel('temp de parcours (ns)');
xlabel('Distance (m)');
title('raw data')
fancy_figure()

% ------------------------------------------
%   cook raw data. remove wow and time shift
% ------------------------------------------

% remove wow
%
d = dewow(d,dt/1000,f0);
% d = dewow_filt(d,dt/1000,f0);

% shift time up
%
delt0 = 15.4; % 1
dels = round(delt0/dt);
d = t0corr(d,dels);

figure;
amp = max(max(d));
imagesc(r,t,d,[-sc*amp sc*amp])
colormap(clmap)
ylabel('temp de parcours (ns)');
xlabel('Distance (m)');
title('cooked data')
fancy_figure()

% ------------------------------------------
%   fourier
% ------------------------------------------

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

subplot(1,2,1)
amp = max(max(d_pow));
imagesc(f,r,d_pow.',[-sc*amp sc*amp])
colormap(clmap)
ylabel('r [m]')
xlabel('f [GHz]')
title('data in frequency')
fancy_figure()

subplot(1,2,2)
hold on
for i=1:nr
plot(f,d_pow(:,i),'.-')
end
hold off
axis tight
xlabel('f [GHz]')
ylabel('d power')
title('data in frequency')
fancy_figure()

% ------------------------------------------
%   cook raw data. for "interpretation"
% ------------------------------------------

% apply power gain
%
gpow = 2.5;
d = gain_rt(d,gpow);
d_gain = d;

% % apply bandpass filter
% %
% f1 = 0.05;      % set low cut
% f2 = 0.1;      % set low shoulder
% f3 = 0.3;     % set hi shoulder
% f4 = 0.38;     % set hi cut
% d = bpfilter_bradford(d,dt,f1,f2,f3,f4);

f_low = 0.05;
f_high = 0.3;
nbutter = 6;
d = butter_butter(d,dt,f_low,f_high,nbutter);

% ------------------------------------------
%   see
% ------------------------------------------

figure;
amp = max(max(d_gain));
imagesc(r,t,d_gain,[-sc*amp sc*amp])
colormap(clmap)
ylabel('temp de parcours (ns)');
xlabel('Distance (m)');
title('overcooked data. gain')
fancy_figure()

figure;
amp = max(max(d));
imagesc(r,t,d,[-sc*amp sc*amp])
colormap(clmap)
ylabel('temp de parcours (ns)');
xlabel('Distance (m)');
title('overcooked data. gain + bandpass')
fancy_figure()
