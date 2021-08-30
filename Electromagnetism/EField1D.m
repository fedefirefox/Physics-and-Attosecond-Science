function [E,Ef,f] = EField1D(FWHM,f0,GDD,TOD,CEP,timeaxis,dt)

% **** Electric Field and Vector Potential Construction from their spectra (Half spectrum!) ****
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

%%%%%%%%%%%%%%%%%%%%
%   CONVERSIONS    %
%%%%%%%%%%%%%%%%%%%%
fmax=1/(2*dt);                    % Maximum Frequency Detectable
F=1/((length(timeaxis)-1)*dt);    % Frequency Resolution (Augmented)
f=0:F:2*fmax;   
sigmat=1/(2*sqrt(2*log(2)))*FWHM; % sigmat in amplitude
sigmaf=1/(sigmat*2*pi);


% Spectrum (Standard sin spectrum) 
Ef=(exp(-(f-f0).^2/(2*sigmaf^2))*exp(-1i*CEP).*exp(-0.5*1i*(2*pi*(f-f0)).^2*GDD).*exp(-1/6*1i*TOD*(2*pi*(f-f0)).^3)); % Time Signal
E=real(fftshift(ifft((Ef))));
% Normalization
norm=max(E);
E=E/norm; 




end

