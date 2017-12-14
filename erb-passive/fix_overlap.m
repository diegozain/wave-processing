function [d,t,i_overlap,i_gap] = fix_overlap(d,t,it_t1,it_t2,i_overlap,i_gap,dt)

% ----------------------------------------------
%   remove & fill overlap
% ----------------------------------------------

% interpolate gap in overlap,
%
% d( t_ ) = interp1(t,d(t),t_)
%
i_ol2gap = [it_t2 , i_overlap(1) + it_t1];
t_ol2gap = t(i_ol2gap);
d_ol2gap = d(i_ol2gap);
t_ol2gap_ = t_ol2gap(1):dt:t_ol2gap(end);
d_ol2gap_ = interp1(t_ol2gap,d_ol2gap,t_ol2gap_);
d_ol2gap_(1) = []; d_ol2gap_(end) = [];
t_ol2gap_(1) = []; t_ol2gap_(end) = [];

% dedekind cut on overlap-gap
%
t_1 = t(1:i_ol2gap(1));
t_2 = t(i_ol2gap(2):end);
d_1 = d(1:i_ol2gap(1));
d_2 = d(i_ol2gap(2):end);

% fill overlap-gap
%
t = [t_1; t_ol2gap_.'; t_2];
d = [d_1; d_ol2gap_.'; d_2];

% ----------------------------------------------
%   fix shifted gaps
% ----------------------------------------------

% # of elements of t (and d) that were removed,
%
n_ol2gap = numel(it_t2+1 : i_overlap(1) + it_t1-1);
% # of elements that *would* have been removed if 
% t1 and t2 were discretized the same,
% numel( t(it_t2+1) :dt: t(i_overlap(1) + it_t1-1) );

% # of elements of t (and d) that were added,
%
m_ol2gap = numel(t_ol2gap_);

i_gap = i_gap - n_ol2gap + m_ol2gap;
i_overlap = i_overlap - n_ol2gap + m_ol2gap;

end