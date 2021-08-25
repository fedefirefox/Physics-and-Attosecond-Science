% Self Phase Modulation Software 
% For didactics and Simulations
%
%
%
%
%
% Politecnico di Milano,
% Federico Vismarra, Attosecond Research Center 2021
% federico.vismarra@polimi.it
%
%
% For further information ask for the guide document.
close all
clc
clear all

%% Constant
c0= 2.99792458e8;         %(m/s)
eps0= 8.8541878128e-12;   
%% FIELD DEFINITION.
FWHM_IR=30;           % IR Duration  (fs)                                                              %
lambda_IR=800;        % Wavelength IR in (nm)                                                          %
GDD=0;                % Input Second Order of Dispersion  (fs^2)                                             %
TOD=0;                % Input Third Order of Dispersion   (fs^3)                                             %
CEP=0;                % (0-2*pi) shift with respect the center                                         %
Ipeak=3e12;           % Peak Intensity in W/cm^2   
f0=c0/(lambda_IR*10^(-9))*10^(-15); % 1/fs 
%% Medium paramters: 
% % % Propagation in air:
% % b2 = 2.1e1; % fs2/m 
% % b3 = 9.5; % fs3/m 
% % n2 = 3e-23;% m2/W
% % n0 = 1.00027505;
% % n2 = (e0*c*n2)/(2*n0);%m^4/V^2
% % alfa = 0.78e-5; %m^-1
% % gamma =(n2*2*pi*f0*10^15)/c;
% 
% % Propagation in fused silica:
 b2 = 35.2e3; % fs2/m     
 b3 = 27.5e3*0; % fs3/m 
 n2 = 2.6e-20;% m2/W
 n0 = 1.45332; 
 n2 = (eps0*c0*n2)/(2*n0);%m^4/V^2
 % alfa = 0.5e-3; %m^-1
 gamma =(n2*2*pi*f0*10^15)/c0;

%% TIME AND FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.002; %Time Resolution (fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeaxis=-5500+dt/2:dt:5500-dt/2;     %Time axis (fs)  ALWAYS EVEN

E=zeros(length(timeaxis),1);
E0=sqrt(2*Ipeak*10^4/(eps0*c0));
%%%%%%%%%%%%%%%%%%%%
%   CONVERSIONS    %
%%%%%%%%%%%%%%%%%%%%
fmax=1/(2*dt);                    % Maximum Frequency Detectable
F=1/((length(timeaxis)-1)*dt);    % Frequency Resolution (Augmented)
f=-fmax:F:fmax;   
sigmat=1/(2*sqrt(2*log(2)))*FWHM_IR; % sigmat in amplitude
sigmaf=1/(sigmat*2*pi);

% Spectrum (Standard sin spectrum) 
Ef=exp(-(f-f0).^2/(2*sigmaf^2)).*exp(-0.5*1i*(2*pi*(f-f0)).^2*GDD).*exp(-1/6*1i*TOD*(2*pi*(f-f0)).^3);
% Time Signal
E=real(fftshift(ifft(ifftshift(Ef))));
% Normalization
norm=max(E);
E=E0*E/norm;


figure(1)
plot(timeaxis,E);
title('input Field')
Aenvelope=E0/norm*exp(-f.^2/(2*sigmaf^2)).*exp(-0.5*1i*(2*pi*f).^2*GDD).*exp(-1/6*1i*TOD*(2*pi*f).^3);
hold on
A=fftshift(ifft(ifftshift(Aenvelope)));
plot(timeaxis,abs(A));
axis([-20,20,-max(abs(A)),max(abs(A))])
ylabel('Electric Field Amplitude (V/m)')
xlabel('Time (fs)')

%%  SPLIT STEP METHOD
figure(2)
A0=ifftshift(fft(fftshift(A)));
L=1e-3; %m
dz=1e-5; %m
A0=A;
 figure(2)
 plot(timeaxis,abs(A0));  
 axis([-20,20,-max(abs(A)),max(abs(A))])
for z=0:dz:L
    figure(3)
   
  
    A1=ifft(ifftshift(exp(-1i*dz/2*(0.5*b2*(2*pi*f).^2+1/6*b3*(2*pi*f).^3)).*fftshift(fft(A0))));    
    
 
  
   A2=exp(-1i*dz*gamma*abs(A1).^2).*A1;
 
    A3=ifft(ifftshift(exp(-1i*dz/2*(0.5*b2*(2*pi*f).^2+1/6*b3*(2*pi*f).^3)).*fftshift(fft(A2))));    
    
    A0=A3;   
    
    % PLOTTING
    subplot(2,2,[1,2])
    plot(timeaxis,abs(A0),'k'); 
    hold on
    plot(timeaxis,-abs(A0),'k'); 
    xlabel('Time (fs)')
    ylabel('Amplitude (V/m)')
    plot(timeaxis,real(A0.*exp(1i*2*pi*f0*timeaxis)),'r')
    hold off

    axis([-50,50,-1.1*max(abs(A)),1.1*max(abs(A))])
    
    subplot(2,2,3)
    AMPSPECTRUM=abs(fftshift(fft(A0)));
    
    plot(f+f0,AMPSPECTRUM/max(AMPSPECTRUM))
    xlabel('Frequency (PHz)')
    ylabel('Arb. units')
    axis([-10*sigmaf+f0,10*sigmaf+f0,0,1])
    
    subplot(2,2,4)
    
    PHASETime=smooth(unwrap(angle(A0)),20);
    N=round(length(PHASETime)/2);
    PHASETime=PHASETime-PHASETime(N);
    instantF=diff(PHASETime)/dt;

    plot(timeaxis(2:end),instantF)
    R=min(find(timeaxis>50));
    L=min(find(timeaxis>-50));
    axis([-50,50,min(instantF(L:R)),max(instantF(L:R))])
     xlabel('Time (fs)')
     ylabel('Relative frequency variation (PHz)')
z   
pause(0.01)
end

