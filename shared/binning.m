function ix_r = binning(x,r)
	
%
% bin numbers x onto receivers r.
%
% does this by returning indicies on x 
% that map x into r:
%
% binned x on r = x( ix_r ), i.e.
%
% 				x(ix_r) ~ r
%
% IMPOrTANT: assumes x and r are sorted.

nx  = length(x);
nr = length(r);

ix_r = ones(nr,1);
dista_ = Inf;
j = 1;
j_ = j;

for i=1:nr
	bol = 1;
	while bol == 1
		dista = abs( r(i)-x(j) );
		if dista_ <= dista
			bol = 0;			% go to next i -> exit while
			ix_r(i) = j_; 		% record index for x
			dista_ = Inf;		% reset distance
			j_ = j;				% record current position
		elseif dista_ > dista
			dista_ = dista;
			j_ = j;
			if j<nx
				j = j+1;
			end
		end
	end
end

end