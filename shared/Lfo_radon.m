function Lfo = Lfo_radon(x,sx,fo)
	
	% discretizes
	%
	% d(x,fo) = sum( exp(...) * p(slow,fo) )_{slow_min,slow_max}
	% 
	% as 
	%
	% d = Lfo * p, Lfo = Lfo(x,slow,fo)
	% to w space
	%
	wo = 2*pi*fo;
	
	% build k space using outer product
	%
	t_shift = 1i*wo * (x' * sx);
	
	% build Lfo
	%
	Lfo = exp(t_shift);

end