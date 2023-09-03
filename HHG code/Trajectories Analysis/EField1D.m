function [A,E] = EField1D(FWHM,f0,GDD,TOD,CEP,timeaxis,dt)

% **** Electric Field and Vector Potential Construction from their spectra ****
%
% FWHM: is the full width at half maximum wanted by the user (au)
% omega: is the central frequency(rad/au)
% GVD: Second Order of Dispersion (au^2)
% TOD: Third Order of Dispersion  (au^2)
% timeaxis: temporal axis wanted by the user (au)
% dt: temporal resolution wanted by the user (au)
%
%
% Author: F.Vismarra (Politecnico of Milano, Attosecond Research Center)
%%


A=zeros(length(timeaxis),1);
E=zeros(length(timeaxis),1);

%%%%%%%%%%%%%%%%%%%%
%   CONVERSIONS    %
%%%%%%%%%%%%%%%%%%%%
fmax=1/(2*dt);                    % Maximum Frequency Detectable
F=1/((length(timeaxis)-1)*dt);    % Frequency Resolution (Augmented)
f=-fmax:F:fmax;   
sigmat=1/(2*sqrt(2*log(2)))*FWHM; % sigmat in amplitude
sigmaf=1/(sigmat*2*pi);


% Spectrum (Standard sin spectrum) 
Af=0.5*(exp(-(f-f0).^2/(2*sigmaf^2))*exp(1i*pi/2+1i*CEP).*exp(-0.5*1i*(2*pi*(f-f0)).^2*GDD).*exp(-1/6*1i*TOD*(2*pi*(f-f0)).^3))+0.5*(exp(-(f+f0).^2/(2*sigmaf^2))*exp(-1i*pi/2-1i*CEP).*exp(0.5*1i*(2*pi*(f+f0)).^2*GDD).*exp(1/6*1i*TOD*(2*pi*(f+f0)).^3));
% Time Signal
A=real(fftshift(ifft(ifftshift(Af))));

E(1:end-1)=-diff(A)/dt;

% Normalization
norm=max(E);

E=E/norm;
A=A/norm; 


end

