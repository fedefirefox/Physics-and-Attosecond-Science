function [A,E] = EField1D_train(FWHM,f0,GDD,TOD,CEP,timeaxis,dt,freq_IR,Npulses,odd,E_IR)

% **** Electric Field and Vector Potential Construction from their spectra for a train pulses****
%
% FWHM: is the full width at half maximum wanted by the user (au)
% omega: is the central frequency(rad/au)
% GVD: Second Order of Dispersion (au^2)
% TOD: Third Order of Dispersion  (au^2)
% timeaxis: temporal axis wanted by the user (au)
% dt: temporal resolution wanted by the user (au)
% freq_IR: driving IR filed frequency
% Npulses: number of pulses for the train
%
% Authors: F.Vismarra, Yingxuan Wu (Politecnico of Milano, Attosecond Research Center)
%%


A=zeros(size(timeaxis));
E=zeros(size(timeaxis));
[E_IRupper,E_IRlower] = envelope(E_IR);
%%
%%%%%%%%%%%%%%%%%%%%
%   CONVERSIONS    %
%%%%%%%%%%%%%%%%%%%%
fmax=1/(2*dt);          % Maximum Frequency Detectable in a.u.
F=1./((length(timeaxis)-1)*dt);   % Frequency Resolution (Augmented)
f=-fmax:F:fmax;   
sigmat=1/(2*sqrt(2*log(2)))*FWHM;
sigmaf=1/(sigmat*2*pi);

% pulse every half period of laser
tau_l= 1/freq_IR;

% Spectrum (Standard sin spectrum)
Af=0.5*(exp(-(f-f0).^2/(2*sigmaf^2)).*exp(1i*pi/2+1i*CEP).*exp(-0.5*1i*(2*pi*(f-f0)).^2*GDD).*exp(-1/6*1i*TOD*(2*pi*(f-f0)).^3)).*exp(2*pi*1i*f*tau_l/4*odd)...
+0.5*(exp(-(f+f0).^2/(2*sigmaf^2))*exp(-1i*pi/2-1i*CEP).*exp(0.5*1i*(2*pi*(f+f0)).^2*GDD).*exp(1/6*1i*TOD*(2*pi*(f+f0)).^3)).*exp(2*pi*1i*f*tau_l/4*odd);
%%

% Even Shape
if odd==0
A=real(fftshift(ifft(ifftshift(Af))));

% Multiplicative factor
n=length(timeaxis);

for j=1:Npulses
Factor(j)=E_IRupper(find(timeaxis<(tau_l/2*j+tau_l/4*odd*j), 1, 'last' ))/max(E_IRupper);    
end


for ij=1:Npulses %add pulse left and right

% time shift using fourier transform properties 
%F{x(t+t0)}= X(i*w)*exp(i*w*t0)
Af_temp= Af.*exp(2*pi*1i*f*ij*tau_l/2)+Af.*exp(-2*pi*1i*f*ij*tau_l/2);

% Time Signal
% as the sign switches for two consecutive pulses *1i^(2n)
A_temp=real(fftshift(ifft(ifftshift(Af_temp),'symmetric')))*1i^(2*ij)*Factor(ij);
A=A+A_temp;
end
toc
E(1:end-1)=-diff(A)/dt;

% Normalization
norm=max(E);

E=E/norm;
A=A/norm; 


figure
subplot(211)
plot(E)
subplot(212)
plot(A)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Odd Shape
% NOTE: THE SHAPING FOR THE ODD MUST BE IMPROVED GAUSSIAN PROFILING OF THE WHOLE PULSE IS NOT
% THE BEST WAY TO PROCEED

if odd==1
A=real(fftshift(ifft(ifftshift(Af),'symmetric')));

% Multiplicative factor
n=length(timeaxis);



for ij=1:Npulses %add pulse left and right

% time shift using fourier transform properties 
%F{x(t+t0)}= X(i*w)*exp(i*w*t0)
Af_temp= Af.*exp(2*pi*1i*f*ij*tau_l/2)+Af.*exp(-2*pi*1i*f*ij*tau_l/2);

% Time Signal
% as the sign switches for two consecutive pulses *1i^(2n)
A_temp=real(fftshift(ifft(ifftshift(Af_temp),'symmetric')))*1i^(2*ij);
A=A+A_temp;
end

E(1:end-1)=-diff(A)/dt;

% Normalization
norm=max(E);

E=E/norm;
A=A/norm; 


end









end

