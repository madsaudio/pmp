function xc=analytic(xr);
% analytic.m
%
% By Mads G. Christensen (christensen@ieee.org)
%
% Calculates the downsampled analytic signal Xc = Xr + i*Xi using the
% FFT, where Xi is the Hilbert transform of the real input signal Xr.
%
% Syntax: Xc=analytic(Xr)
%

X=fft(xr);
X(1)=X(1)./2;
X(end/2+1)=X(end/2+1)/2;
X(end/2+2:end)=0;
xc=2*ifft(X);
xc=xc(1:2:end-1);