clear all;
close all;
clc;

path='';
filename='abba.wav';
[input,setup.fs]=audioread([path filename]);
input=input(:,1);
input=input*2^15;
output=zeros(length(input),1);

Q=1;
input=resample(input,1,Q);
setup.fs=setup.fs/Q;

setup.F=2^16;                       % FFT size
setup.N=2*round(0.030*setup.fs/2);    % Segment length
setup.w=hanning(setup.N);             % Analysis/Synthesis window
setup.L=30;                           % Number of sinusoids
setup.I=10;

psy=psy_init(setup.fs,setup.F);

n=0:setup.N-1;
N=setup.N;
F=setup.F;
w=setup.w;
pos=1:setup.N;

wb=waitbar(0,'Processing...');
while max(pos) <= length(input),
    x=(input(pos));
    if var(x)>100,
        [psy.msk] = psy_comp(x,w,psy);
        wei=1./(psy.msk+eps);      
        [a_pmp,f_pmp,z_pmp,d_pmp]=analysis_pmp(x,wei,setup);
        xhat=real(z_pmp);
        output(pos)=output(pos)+xhat;
        
        
        figure(1);clf,plot(20*log10(abs(fft(x,F))));hold on
        plot(20*log10(psy.msk));
        pause
        
        
    end    
    pos=pos+(N/2);
    waitbar(pos(end)/length(input),wb);
end;
close(wb);
fileout=['output_' filename];
audiowrite(fileout,output*2^-15,setup.fs);


