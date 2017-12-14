function [i_overlap,i_gap] = overlaps_gaps(X,I)

% get all overlaps and gaps from miniSEED 
% data read with 
% rdmseed.m,
%
% [X,I] = rdmseed( filename );


% ------------------------------------------
% overlaps
% ------------------------------------------

i_overlap_block = I.OverlapBlockIndex;
n_overlap_block = numel( i_overlap_block );
i_overlap = zeros(n_overlap_block,1);

% i_overlap will have indicies of d and t 
% where the overlaps begin, 
% 
% i.e. indicies of j'th overlap are
% 
% i_overlap(j) --- i_overlap(j) + 1
% 
% j'th overlap is
% 
% d(i_overlap(j)) --- d(i_overlap(j) + 1)
% t(i_overlap(j)) --- t(i_overlap(j) + 1)

for j=1:n_overlap_block
id_overlap = 0;
for i=1:i_overlap_block(j)-1
id_overlap = id_overlap + X(i).NumberSamples;
end
i_overlap(j) = id_overlap;
end

% ------------------------------------------
% gaps
% ------------------------------------------

i_gap_block = I.GapBlockIndex;
n_gap_block = numel( i_gap_block );
i_gap = zeros(n_gap_block,1);

% same idea as with overlaps.

for j=1:n_gap_block
id_gap = 0;
for i=1:i_gap_block(j)-1
id_gap = id_gap + X(i).NumberSamples;
end
i_gap(j) = id_gap;
end

end