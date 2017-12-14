function [d_fk,f,k,df,dk] = fk_rt(d,dt,dx)

% output is (nf x nk)

[nt,nx] = size(d);
df = 1/dt/nt;
dk = 1/dx/nx;

% fft
%
d_fk = fft(fft(d).').';

f = (-nt/2:nt/2-1) * df;
k = (-nx/2:nx/2-1) * dk;

% fft shift
%
d_fk = fftshift(d_fk,1);
d_fk = fftshift(d_fk,2);

% get rid of negative part
%
d_fk = d_fk( 1:ceil(nt/2)-1 , : );
d_fk = d_fk( : , ceil(nx/2)+1:nx-1);

f = f( ceil(nt/2)+1:nt-1 );
k = k( ceil(nx/2)+1:nx-1 );

d_fk = flip(d_fk);

end
