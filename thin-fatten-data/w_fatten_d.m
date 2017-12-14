function d_fat = w_fatten_d(d_thin,x,r,t_)

%
% old - (received)
%
[R,T_] = meshgrid(r,t_);
%
% new - (for synthetics)
%
[X,T] = meshgrid(x,t_);

d_fat = interp2(R,T_,d_thin,X,T);

end