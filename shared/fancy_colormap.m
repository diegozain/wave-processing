function [clmap,c_axis,amp] = fancy_colormap(a)

%
% a is a matrix to be plotted using imagesc 
%

% for plotting
%
sc=0.3;
clo=[0 0 1];cmid=[.99 .99 .99];chi=[.5 .12 0];
clmap=mkclrmap(clo,cmid,chi);

amp = max(a(:));
c_axis = [-sc*amp sc*amp];

end