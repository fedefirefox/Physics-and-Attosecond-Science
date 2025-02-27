%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WELCOME TO 1d RECOVERY TESTING ENVIRONMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v.2.0
%Ricostruzione di fase monodimensionale con fft, ricostruisco gaussiana con ampiezze da 0.01 a 10__ LA PRECISIONE VARIA CON ;
%OTTIMIZZATA PER PROFILO DI FASE GAUSSIANO 0.1-2 di FWHM range

% QUI DEFINISCI I PARAMETRI DEL TUO OGGETTO MONODIMENSIONALE SIMULATO
N=1485;
dx=0.0044; %mm one pixel is 4.4um
x=-N/2*dx:dx:N/2*dx-dx;
phase=0.1*(exp(-100*x.^2)); 
figure(5)
plot(x,phase);
E0=1;
vx=7; 
extrashift=0.2;
Itot=(2*E0+E0*cos(2*pi*vx*x+phase+extrashift));

%Itot=fringesLighter(500,:);
%N=length(Itot);
%x=-N/2*dx:dx:N/2*dx-dx;
%% Reconstruction via fourierTransform.
dFx=1/(N*dx); %Minima(x,y) spatial frequencies detectable
Fsx=1/dx; % almost max frequency detectable 
fx=(-Fsx/2:dFx:Fsx/2-dFx) %As Usual for the FFT 
IF=fftshift(fft(Itot));
figure(1)
plot(fx,abs(IF));
xlabel('XSpacialFreq(1/mm)')
ylabel('Amplitude')
title("AMPLITUDE ");
IFm=abs(IF)
IFm(1:round(N/2)+5)=0;
sigma=3;
figure(2)
plot(fx,IFm)
xlabel('XSpacialFreq(1/mm)')
ylabel('Amplitude')
title("AMPLITUDE Modified ");

C = xcorr(exp(-(fx.^2/(2*sigma^2))).^40,IFm);
[fxnot]=find(C==max(C))
fxPEAK=-(fxnot-N)*dFx+dFx/2
sigma2=4.5;
IF2=IF.*exp(-((fx-fxPEAK).^2/(2*sigma2^2))).^20;

figure(3)
plot(fx,abs(IF2))
xlabel('XSpacialFreq(1/mm)')
ylabel('Amplitude')
title("AMPLITUDE Filtered ");

IT=ifft(ifftshift(IF2));

ReconstructedPhasis=unwrap(angle(IT.*exp(-1i*2*pi*vx*x)));
% BUTTA VIA I BORDI
ReconstructedPhasis(1:round(N/2-2.5/dx))=0;
ReconstructedPhasis(round(N/2+2.5/dx):end)=0;

% With a gaussian profile many bad shit occurs with these trick i can

% restore the flatness 1) fight gaussian profile distortion

for i=1:5
m=(ReconstructedPhasis(round(N/2+2.5/dx)-i)-ReconstructedPhasis(round(N/2-2.5/dx)+i))/(x(round(N/2+2.5/dx)-i)-x(round(N/2-2.5/dx)+i));
ReconstructedPhasis=ReconstructedPhasis-m*x;
end
%restore the right height
for i=1:5
h=(ReconstructedPhasis(round(N/2+2.5/dx)-i)+ReconstructedPhasis(round(N/2-2.5/dx)+i))/2
ReconstructedPhasis=ReconstructedPhasis-h;
end
ReconstructedPhasis(1:round(N/2-2.5/dx))=0;
ReconstructedPhasis(round(N/2+2.5/dx):end)=0;


figure(4)
plot(x,ReconstructedPhasis);
xlabel('X(mm)');
ylabel('PhaseShift');
title("Reconstructed Phase Profile");