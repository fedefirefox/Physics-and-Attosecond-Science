%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Welcome in PhaseRecovery v2.0 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v.2.0 24/10/2019
% AIM: ACHIEVE GOOD RESULT with vx=10


% This algorithm reconstruct the hidden phase of a quasi-monochromatic 1D object
% Itot=E0(1+cos(2*pi*vx*x+phase(x)));

% The algorithm works well if: vx=10, phase=Aexp(-Bx^2) A=0.01-10 and B=0.1-100 for "Fourier reason"

%You can try SimulationPhaseRECO1D for testing

%You can try FinaleSimulationFringes for 2D interference fringes recovery
function PHASE=PhaseRecovery(Itot,x,dx,N)

%Extension of Itot to increase frequency resolution?

%Definition of Fourier Parameters:
dFx=1/(N*dx); %Minimum(x spatial frequencies detectable
Fsx=1/dx; % almost max frequency detectable 
fx=(-Fsx/2:dFx:Fsx/2-dFx); %As Usual for the FT 

IF=fftshift(fft(Itot)); % Direct Fourier %If it does not work extend the number of points in Itot!!!

% Automatic detection of the peak at positive frequency
IFm=abs(IF);
IFm(1:round(N/2)+20)=0;
sigma=3;
C = xcorr(exp(-(fx.^2/(2*sigma^2))).^40,IFm); %superGaussian cross-correlation
[fxnot]=find(C==max(C));
fxPEAK=-(fxnot-N)*dFx+dFx/2;


%Recomputing everthing with new x-reshifted 


sigma2=4.5;
%Filtering only of the right part for phaseRecovery
IF2=IF.*exp(-((fx-fxPEAK).^2/(2*sigma2^2))).^6; %SuperGaussian filtering

IT=ifft(ifftshift(IF2));

ReconstructedPhase=unwrap(angle(IT.*exp(-1i*2*pi*fxPEAK*x))); % MIGHT Here is essentiatial the centering of the function
%Remove the boundary that are ugly

%2 is an arbitrary decision
death=2.5;
ReconstructedPhase(1:round(N/2-death/dx))=0;
ReconstructedPhase(round(N/2+death/dx):end)=0;
PHASE=ReconstructedPhase;

end

%"Fourier Reason" See https:..... fedepaper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Federico Vismarra 19/10/2019 POLIMI ATTOSECONDLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%