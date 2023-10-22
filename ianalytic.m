function xr=ianalytic(xc);
% ianalytic.m
%
% By Mads G. Christensen (christensen@ieee.org)
%
%
% Syntax: Xr=ianalytic(Xc)
%

if ~(size(xc,2)==1),
  disp('Error')
  return
end

X=fft(xc);
X=[X; zeros(size(X))];
xr=2*real(ifft(X));
