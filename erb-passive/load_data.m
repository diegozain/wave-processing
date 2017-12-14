
% ---------------------------------------------
% loads miniseed data into matlab
% ---------------------------------------------

% all these files are in ERB_data/XX/etc

% Stations are
% BSM1 (center),
% BSM3 (west),
% BSM7 (east, alumni center station)

BSM1_1 = 'ERB_data/XX/BSM1/XX.BSM1..HHX_MC-PH1_0429_20171201_160000.miniseed';
BSM1_2 = 'ERB_data/XX/BSM1/XX.BSM1..HHX_MC-PH1_0429_20171201_170000.miniseed';
BSM1_3 = 'ERB_data/XX/BSM1/XX.BSM1..HHY_MC-PH1_0429_20171201_160000.miniseed';
BSM1_4 = 'ERB_data/XX/BSM1/XX.BSM1..HHY_MC-PH1_0429_20171201_170000.miniseed';
BSM1_5 = 'ERB_data/XX/BSM1/XX.BSM1..HHZ_MC-PH1_0429_20171201_160000.miniseed';
BSM1_6 = 'ERB_data/XX/BSM1/XX.BSM1..HHZ_MC-PH1_0429_20171201_170000.miniseed';

BSM3_1='ERB_data/XX/BSM3/XX.BSM3..HHX_MC-PH1_0423_20171201_160000.miniseed';
BSM3_2='ERB_data/XX/BSM3/XX.BSM3..HHX_MC-PH1_0423_20171201_170000.miniseed';
BSM3_3='ERB_data/XX/BSM3/XX.BSM3..HHY_MC-PH1_0423_20171201_160000.miniseed';
BSM3_4='ERB_data/XX/BSM3/XX.BSM3..HHY_MC-PH1_0423_20171201_170000.miniseed';
BSM3_5='ERB_data/XX/BSM3/XX.BSM3..HHZ_MC-PH1_0423_20171201_160000.miniseed';
BSM3_6='ERB_data/XX/BSM3/XX.BSM3..HHZ_MC-PH1_0423_20171201_170000.miniseed';

BSM7_1='ERB_data/XX/BSM7/XX.BSM7..HHX_MC-PH1_0426_20171201_160000.miniseed';
BSM7_2='ERB_data/XX/BSM7/XX.BSM7..HHX_MC-PH1_0426_20171201_170000.miniseed';
BSM7_3='ERB_data/XX/BSM7/XX.BSM7..HHY_MC-PH1_0426_20171201_160000.miniseed';
% BSM7_4='ERB_data/XX/BSM7/XX.BSM7..HHY_MC-PH1_0426_20171201_170000.miniseed.part'
BSM7_5='ERB_data/XX/BSM7/XX.BSM7..HHZ_MC-PH1_0426_20171201_160000.miniseed';
BSM7_6='ERB_data/XX/BSM7/XX.BSM7..HHZ_MC-PH1_0426_20171201_170000.miniseed';

% ---------------------------------------------
% get usefull info
% ---------------------------------------------
% for more info do:
% 
% fieldnames(X)
% fieldnames(I)

% ------------------------------------------------
% center station
% 
% BSM1
% ------------------------------------------------
clear d t X I

[X,I] = rdmseed(BSM1_5);

% -------------------------------------------------------
%
%   example
%
% -------------------------------------------------------

t = cat(1,X.t) * (24*60*60);
d = cat(1,X.d);

fs = X.SampleRate; % [Hz]
dt=1/fs; % [s] this is in seconds because of *(24*60*60)

% ----------------------------------------------
%   get overlap and gap intervals
% ----------------------------------------------

% for i'th overlap (similar for gap)
% the two points determining it are,
%
% [ t( i_overlap(i) ),d( i_overlap(i) ) ] -- 
%  -- [ t( i_overlap(i)+1 ),d( i_overlap(i)+1 ) ]
%
[i_overlap,i_gap] = overlaps_gaps(X,I);

% ----------------------------------------------
%   get gap-interval in overlap 
% ----------------------------------------------

% for the i'th overlap
% the two points that will determine the gap are,
% 
% [ t1(it_t2),d(it_t2) ] -- [ t2(it_t1),d(i_overlap(1) + it_t1) ]
%
%
% on t this is, 
%
% [ t(it_t2),d(it_t2) ] -- 
% -- [ t(i_overlap(1) + it_t1),d(i_overlap(1) + it_t1) ]
%
t1 = t(1:i_overlap(1));
t2 = t(i_overlap(1)+1:end);
[it_t2,it_t1] = get_overlap(t1,t2);

% -------------------------------------------------------
%
%   see raw data
%
% -------------------------------------------------------

figure;
hold on
% ------------------
% data
% ------------------
plot(t,d,'.-')
% ------------------
% overlaps
% ------------------
for i=1:numel(i_overlap)
plot(t(i_overlap(i)),d(i_overlap(i)),'.g','Markersize',30)
plot(t(i_overlap(i)+1),d(i_overlap(i)+1),'.r','Markersize',30)
end
% ------------------
% gaps
% ------------------
for i=1:numel(i_gap)
plot(t(i_gap(i)),d(i_gap(i)),'sg','Markersize',10,...
'MarkerFaceColor','g')
plot(t(i_gap(i)+1),d(i_gap(i)+1),'sr','Markersize',10,...
'MarkerFaceColor','r')
end
% ------------------
% overlap interval
% ------------------
plot(t1( it_t2 ),d(it_t2),'^g','Markersize',10,...
'MarkerFaceColor','g')
plot(t2( it_t1 ),d(it_t1+i_overlap(1)),'^r','Markersize',10,...
'MarkerFaceColor','r')
hold off
% datetick('x','keeplimits')
grid on
grid minor
axis tight
xlabel('$t$ [s]')
ylabel('counts')
title('$d(t,\mbox{center})$ rawest')
fancy_figure()

% -------------------------------------------------------
%
%   FIX
%
% -------------------------------------------------------

% ----------------------------------------------
%   fix overlap
% ----------------------------------------------

[d,t,i_overlap,i_gap] = fix_overlap(d,t,it_t1,it_t2,i_overlap,i_gap,dt);

% ----------------------------------------------
%   fix gap
% ----------------------------------------------

[d,t] = fix_gap(d,t,i_gap,dt);

figure;
% ------------------
% data
% ------------------
plot(t,d,'.-')
grid on
grid minor
axis tight
xlabel('$t$ [s]')
ylabel('counts')
title('$d(t,\mbox{center})$ fixed but raw')
fancy_figure()

% -------------------------------------------------------
%
%   end example
%
% -------------------------------------------------------

[X,I] = rdmseed(BSM1_6);
[d_,t_,~] = fix_trace(X,I);
t = cat(1,t,t_);
d = cat(1,d,d_);

dz_center = d;
t_center = t;

% figure;
% plot(t,d,'.-')
% grid on
% grid minor
% axis tight
% xlabel('time [s]','FontSize',20)
% ylabel('counts','FontSize',20)
% title('raw data. center station','FontSize',20)

% ------------------------------------------------
% west station
% 
% BSM3
% ------------------------------------------------
clear d t X I

[X,I] = rdmseed(BSM3_5);
[d,t,dt] = fix_trace(X,I);

[X,I] = rdmseed(BSM3_6);
[d_,t_,~] = fix_trace(X,I);

t = cat(1,t,t_);
d = cat(1,d,d_);

dz_west = d;
t_west = t;

% figure;
% plot(t,d,'.-')
% grid on
% grid minor
% axis tight
% xlabel('time [s]','FontSize',20)
% ylabel('counts','FontSize',20)
% title('raw data. west station','FontSize',20)

% ------------------------------------------------
% east station
% 
% BSM7
% ------------------------------------------------
clear d t X I

[X,I] = rdmseed(BSM7_5);
[d,t,dt] = fix_trace(X,I);

[X,I] = rdmseed(BSM7_6);
[d_,t_,~] = fix_trace(X,I);

t = cat(1,t,t_);
d = cat(1,d,d_);

dz_east = d;
t_east = t;

% figure;
% plot(t,d,'.-')
% grid on
% grid minor
% axis tight
% xlabel('time [s]','FontSize',20)
% ylabel('counts','FontSize',20)
% title('raw data. east station','FontSize',20)


% ----
% all stations
% ----

figure;
hold on
plot(t_west,dz_west,'m-')
plot(t_center,dz_center,'-')
plot(t_east,dz_east,'g-')    
hold off
axis tight
xlabel('$t$ [s]')
ylabel('counts')  
legend({'west','center','east'},'Location','best')
title('$d(t,r)$ z channel. erb passive')
fancy_figure()







