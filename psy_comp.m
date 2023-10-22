function msk = psy_comp(x,w,psymod);
%PSY_INIT   Initialize psychoacoustic model
%
% Syntax:
%   msk=psy_comp(x,window,psy);
%
% Input:
%   x        input segment
%   w        window function
%   psy      model struct
%
%   
% Output:
%   msk      masking curve
%
% Description:
%   This function computes the masking curve of the input segment x for a
%   given window. It is based on the psychoacoustic model of the paper van 
%   de Par, S., Kohlrausch, A., Heusdens, R. et al. A Perceptual Model for 
%   Sinusoidal Audio Coding Based on Spectral Integration. EURASIP J. Adv. 
%   Signal Process. 2005, 317529 (2005). It requies that the struct of the 
%   psychoacoustic model, which can be initialized with psy_init(), be 
%   given as input.
%
% Example:
%   msk=psy_comp(x,hanning(length(x),psy);
%
F=psymod.fftsize;
N=length(x);
psi=1e-32;  % 10*log10(msk)<=10*32
fb=psymod.fb;                             
f=[0:1:0.5*F]'/F*psymod.fs;
X=fft(x.*w,F);
Pxx=(abs(X(1:F/2+1))).^2/F/norm(w)^2;	
tmp=Pxx(1:F/2+1)./psymod.thr(1:F/2+1);	
aux=fb.^2*tmp + psymod.C_a;
d=(fb'.^2*(1./aux))*(psymod.C_s)./psymod.thr;
msk = 1./(d+psi);	
msk(end+1:2*(end-1))=flipud(msk(2:end-1));