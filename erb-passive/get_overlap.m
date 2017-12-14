function [it_t2,it_t1] = get_overlap(t1,t2)
	
%
% get index of overlapped t1 on t2 and viceversa:
%
%					it_t2
%						|
% t1= .--.--.--.--.
%						.--.--.--.--. =t2
%						|
%						1
%
% 		t1( it_t2 ) ~ t2(1)
%
%
%
%								t1(end)
%									|
% t1= .--.--.--.--.
%						.--.--.--.--. =t2
%									|
%								it_t1
%
% 		t2( it_t1 ) ~ t1(end)
%
%
% IMPOrTANT: assumes t1 and t2 are sorted.

nt1  = length(t1);
nt2 = length(t2);

% ------------------------------
%					it_t2
%						|
% t1= .--.--.--.--.
%						.--.--.--.--. =t2
%						|
%						1
%
% 		t1( it_t2 ) ~ t2(1)
% ------------------------------

dista_ = Inf;
j = nt1;
j_ = j;
bol = 1;
while bol == 1
	dista = abs( t2(1)-t1(j) );
	if dista_ <= dista
		bol = 0;
		it_t2 = j_;
	elseif dista_ > dista
		dista_ = dista;
		j_ = j;
		j = j-1;
	end
end

% ------------------------------
%								t1(end)
%									|
% t1= .--.--.--.--.
%						.--.--.--.--. =t2
%									|
%								it_t1
%
% 		t2( it_t1 ) ~ t1(end)
% ------------------------------

dista_ = Inf;
j = 1;
j_ = j;
bol = 1;
while bol == 1
	dista = abs( t2(j)-t1(nt1) );
	if dista_ <= dista
		bol = 0;
		it_t1 = j_;
	elseif dista_ > dista
		dista_ = dista;
		j_ = j;
		j = j+1;
	end
end

end