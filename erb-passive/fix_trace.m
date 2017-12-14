function [d,t,dt] = fix_trace(X,I)

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

if numel(i_overlap) > 0
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
  
  % ----------------------------------------------
  %   fix overlap
  % ----------------------------------------------
  [d,t,i_overlap,i_gap] = fix_overlap(d,t,it_t1,it_t2,i_overlap,i_gap,dt);
end  
if numel(i_gap) > 0
  % ----------------------------------------------
  %   fix gap
  % ----------------------------------------------    
  [d,t] = fix_gap(d,t,i_gap,dt);
end

end