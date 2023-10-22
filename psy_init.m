function psymod=psy_init(fs,F);
%PSY_INIT   Initialize psychoacoustic model
%
% Syntax:
%   psy=psy_init(fs,F);
%
% Input:
%   fs       sampling frequency
%   F        FFT size
%
% Output:
%   psy      struct containing psychocacoustic model
%
% Description:
%   This function initializes the psychoacoustic model of the paper van de
%   Par, S., Kohlrausch, A., Heusdens, R. et al. A Perceptual Model for
%   Sinusoidal Audio Coding Based on Spectral Integration. EURASIP J. Adv.
%   Signal Process. 2005, 317529 (2005). It requires the FFT size and the
%   sampling frequency as input and uses that to initialize the required
%   struct.
%
% Example:
%   psy=psy_init(8192,44100);
%
erbw=@(x)((x.*0.00437+1).*24.7);
inverbrate=@(x)(exp(x.*0.1079)-1)./0.00437;
erbrate=@(x)(log(x.*0.00437+1)./0.1079);
gtfilter=@(x,y,z)(1./(1+((x-y)./(16/(pi*5)*z)).^2).^2);

% Init constants
psi=1e-6;
freq=[0:1:F-1]'/F*fs;
nf=60;

% Save constants
psymod.fs=fs;
psymod.fftsize=F;

% Init gammatone filters
delta=erbrate(fs/2)/nf;
tmp=zeros(nf, F);
for m=1:nf,
    fc=inverbrate(m*delta);
    bw=erbw(fc);
    tmp(m,:)=gtfilter(freq,fc,bw)';
end
h=zeros(nf, F+1);
h(:, 1:F) = tmp(:, 1:F);
h(:, F+1) = tmp(:, 1);
fb = zeros(nf,F/2+1);
fb = h(:,1:F/2+1) + h(:,end:-1:F/2+1);

% Compute centerfreq, bw, absolute threshold
f0v=inverbrate(delta.*(1:nf));
bwv=max(0,erbw(f0v));
% model of threshold in quiet
thr = 3.64*((freq(1:F/2+1)+psi)/1e3).^(-0.8)-6.5*exp(-0.6*(freq(1:F/2+1)/1e3-3.3).^2)+1e-3*(freq(1:F/2+1)/1e3).^4;
% simplified model
%ThresholdQuiet=3.64*(freq(1:F/2+1)/1e3+psi+0.5).^(-0.8)-6.5*exp(-0.6*(freq(1:F/2+1)/1e3-3.3).*(freq(1:F/2+1)/1e3-3.3));
thr = 10.^(thr/10);

% Calibrate model & save
[psymod.C_a,psymod.C_s]=aux_calibration(freq(1:F/2+1),fb,thr);

% Save result in struct
psymod.fb=fb;
psymod.thr=thr;

end

function [C_a, C_s] = aux_calibration(f,fb,thr);
[~, ndx]=min(abs(f-1000));  % find index closest to 1 kHz
f1khz=fix((ndx));
masker = 10^7;	   % masker power (70 dB SPL)
maskee = 10^(5.2); % maskee power (52 dB SPL)
spec = masker/thr(f1khz);
C_tmp = 1;
exc1 = (spec*(fb(:, f1khz).^2)) + C_tmp;
C_s= 1/(sum((fb(:,f1khz)).^2.*(maskee/thr(f1khz))./exc1));  % C_s
C_a = (C_s.*sum(fb(:,f1khz).^2)); % C_a
while abs((C_tmp-C_a)./C_tmp)>1e-6;
    exc1 = (spec*(fb(:, f1khz).^2)) + C_a;                
    C_s= 1/(sum((fb(:,f1khz)).^2.*(maskee/thr(f1khz))./exc1));  % C_s
    C_tmp = C_a;
    C_a = (C_s.*sum(fb(:,f1khz).^2));                     % C_a    
end
end