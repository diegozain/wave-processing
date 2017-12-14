function [d_,f,df] = fourier_rt(d,dt)

% -------------------------------------------------------
% fourier: d(r,t), df -> d(r,f), f
% -------------------------------------------------------

[nt,~] = size(d);
df = 1/dt/nt;

% fft
%
d_ = fft(d,[],1);

% fft shift
%
d_ = fftshift(d_,1);
f = (-nt/2:nt/2-1) * df;

% get rid of negative part
%
d_ = d_( ceil(nt/2)+1:nt-1, : );
f = f( ceil(nt/2)+1:nt-1 );

end
