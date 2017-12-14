function [d_,f,df] = fourier_rt2(d,dt)

% -------------------------------------------------------
% fourier: d(r,t), df -> d(r,f), f
% -------------------------------------------------------

[nt,~] = size(d);

nfft   = 2 * ( 2^nextpow2( nt ) ); 
fs     = 1  / dt;                    
df     = fs / nfft;                  
fny   = fs / 2;                     
f = ( 0 : ( nfft - 1 ) ) .* df; 

% set the correct negative part of frequency array
%
f( f >= fny ) = f( f >= fny ) - ( fny * 2 );
f = f( f > 0);

% fft the data along time axis
%
d_ = fft( d, nfft, 1); 

end