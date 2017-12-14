function d_thin = w_thin_d(d,x,r,t_)

% get indicies of x that map to r, i.e.
%
% x(ix_r) ~ r
%
ix_r = binning(x,r);

% % get indicies of t that map to t_, i.e.
% %
% % t(it_t_) ~ t_
% %
% it_t_ = binning(t,t_);

d_thin = d(:,ix_r);

end