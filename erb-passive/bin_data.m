% --------------------------------------------
% bin selected traces
% --------------------------------------------

% t_west( iwest_east ) ~ t_east
%
iwest_east = binning(t_west,t_east);

t = t_west(iwest_east(1) : iwest_east(end));
dz_west_ = dz_west(iwest_east(1) : iwest_east(end));

% t_east( ieast_west ) ~ t_west
%
ieast_west = binning(t_east,t_west(end));

t_ = t_east(1 : ieast_west);
dz_east_ = dz_east(1 : ieast_west);

figure;
hold on
plot(t_west,dz_west,'m-')
plot(t_east,dz_east,'g-')
plot(t,dz_west_,'b-')
plot(t_,dz_east_,'r-')
hold off
axis tight
xlabel('$t$ [s]')
ylabel('counts')  
legend({'west','east','west-binned','east-binned'},'Location','best')
title('$d(t,r)$ z channel. erb passive')
fancy_figure()

t = t(1):dt:t(end);
nt = numel(t);
d = zeros(nt,2);
d(:,1) = dz_west_; d(:,2) = dz_east_;

clear t_ t_east t_west ieast_west iwest_east dz_east dz_west dz_east_ dz_west_
clear i_gap i_overlap X I d_ it_t1 it_t2 i t1 t2

figure;
hold on
plot(t,d(:,1),'b-')
plot(t,d(:,2),'r-')
hold off
axis tight
xlabel('$t$ [s]')
ylabel('counts')  
legend({'west','east'},'Location','best')
title('$d(t,r)$ z channel. erb passive')
fancy_figure()
