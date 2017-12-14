function [disper_g, f, v, dro_snr] = ftan_2(dro,dt,fmin,fmax,df,vmin,vmax,dv,alph,dsr,isPlot)
%
% USAGE: [fvMatrix, f, v, dro_snr] = ftandro(dro,dt,fmin,fmax,df,vmin,vmax,dv,alph,dsr,isPlot)
%
% Apply frequency-time analysis (FTAN) to a seismic dro. All
% preprocessing (e.g., windowing) should be done prior to calling this
% function.
%
% INPUT:
%   dro = data vector; size( 1, nt )
%   dt    = time sample interval in dro [s]
%   fmax  = maximum frequency [Hz]
%   fmin  = minimum frequency [Hz]
%   df    = frequency sampling interval [Hz]  
%   vmax  = maximum group velocity to scan over [m/s]
%   vmin  = minimum group velocity to scan over [m/s]
%   dv    = velocity sampling interval [m/s]  
%   alph = a dimensionless parameter determining the narrowness of the filtering
%   dsr  = source-receiver dsrance [m]
%   isPlot = plotting flag (=0, no plot; =1 to plot filtered dros)
% OUTPUT:
%   fvMatrix = frequency vs. velocity matrix; size( nv, nf )
%   f   = frequency vector for fvMatrix; size( nf, 1 )
%   v   = velocity vector for fvMatrix; size( 1, nv )
%   dro_snr = frequency dependent signal to noise ratio ( nf, 1 )


% check the shape of dro input
if size(dro,1) ~= 1
    dro = transpose( dro );
end

t = (0:numel(dro)-1).*dt;

nt   = numel( dro ); % update number of samples
nfft  = nt; % number of points in FFT: must be the same as nt
w  = 2 .* pi .* makeFFTvector( dt, nfft); % FFT frequency vecotr (rad)
t = ( 0 : nt-1 ) .* dt; % time vector [s]

% make velocity and frequency vectors
f = fmin : df : fmax;  % [Hz] frequencies to scan over
v = vmin : dv : vmax; % [m/s] velocities to consider

nf = numel(f); % number of elements in frequency vector
nv = numel(v); % number of elements in velocity vector

% apply frequency filter 
dro_filt = zeros( nf, nt );

% make sure to taper the edges to zero before FFT
taper = tukeywin(nt,0.1);
dro = dro .* taper';

for ii = 1 : nf % loop through frequencies
    
    wo =  2 * pi * f(ii); % [rad] central frequency of Guassian

    % make the narrowband filter (Gaussian as in Levshin's work)
    gfl = exp( -alph * ( ( ( w - wo ) ./ wo ).^2 ) ); % equation 6 in Bensen et al. (2007)

    dro_fft = fft( hilbert(dro), nfft) ./ nfft; % FFT of analytic signal

    dro_filt( ii, : ) = ifft( dro_fft .* gfl, nfft); % multiply by narrowband filter and inverse FFT
    % Equation 5 in Bensen et al. (2007)
    
end

disper_g = zeros( nv, nf);
% interpolate from time axis to velocity axis using a spline
for ii = 1 : nf
    disper_g( :, ii ) = spline( t, abs( dro_filt( ii, : ) ), dsr ./ v );
end

tmin = dsr/vmax;
tmax = dsr/vmin;
nmin = round(tmin/dt);
nmax = min([round(tmax/dt) nt]); % make sure nmax is not larger than nt
sigIdx = nmin : nmax;
noiIdx = nmax : nt;

% Compute the SNR of the bandpassed dros
dro_snr = zeros( nf, 1 );
for ii = 1 : nf
    dro_snr( ii ) = max( abs( dro_filt( ii, sigIdx ) ) ) /...
        rms( dro_filt( ii, noiIdx ) );
end

end