function [a,f,z,D]=analysis_pmp(s,det,setup)
w=setup.w(1:2:end);
s=w.*analytic(s);
L=setup.L;
F=setup.F;
N=length(s);
R=fft(s,F);
W=fft(w,F);
v=real(ifft(W));
n=[0:N-1]';
D(1)=sum(det.*abs(R).^2)/F;
z=zeros(N,1);
gg=real(fft(fft(abs(W).^2,F).*ifft(det,F),F))/F;
for l=1:L,
    rg=(fft(ifft(det.*R).*v,F));
    J=abs(rg./sqrt(gg)).^2;
    [tmp,ndx]=max(J);
    f(l)=2*pi*(ndx-1)/F;
    a(l)=rg(ndx)./(gg(ndx));
    psindx=shifti(F,ndx-1);
    R=R-(a(l))*W(psindx);
    z=z+a(l).*w.*exp(j*f(l)*n);
    D(l+1)=sum(det.*abs(R).^2)/F;
end
z=ianalytic(z);

end

function [i]=shifti(N,n);
i=[N-n+1:N 1:N-n];
end
