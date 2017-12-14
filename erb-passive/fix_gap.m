function [d,t] = fix_gap(d,t,i_gap,dt)

% ----------------------------------------------
%   remove & fill gaps
% ----------------------------------------------

% interpolate gap,
%
% d( t_ ) = interp1(t,d(t),t_)
%
for i=1:numel(i_gap)
ii_gap = [i_gap(i) , i_gap(i)+1];
t_gap = t(ii_gap);
d_gap = d(ii_gap);
t_gap_ = t_gap(1):dt:t_gap(end);
d_gap_ = interp1(t_gap,d_gap,t_gap_);
d_gap_(1) = []; d_gap_(end) = [];
t_gap_(1) = []; t_gap_(end) = [];

% dedekind cut on gap
%
t_1 = t(1:ii_gap(1));
t_2 = t(ii_gap(2):end);
d_1 = d(1:ii_gap(1));
d_2 = d(ii_gap(2):end);

% fill gap
%
t = [t_1; t_gap_.'; t_2];
d = [d_1; d_gap_.'; d_2];

% ----------------------------------------------
%   fix shifted gaps
% ----------------------------------------------

% # of elements of t (and d) that were added,
%
m_gap = numel(t_gap_);
i_gap(i+1:end) = i_gap(i+1:end) + m_gap;
end

end