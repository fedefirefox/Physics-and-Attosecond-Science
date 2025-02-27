%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING CODE 2D SIMULATION ENVIRONEMENT               v2.0           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction of a simulated 2D phase profile is implemented 
% The code use the function Phase Recovery and perform a row by row fft
%
% You can use SimulationPhaseREC1D.m to test it 
%
% In the codes are present optimization parameters.
%
% The simulation has been tested for
%
%  Refractive index magnitude                    Dn0=[e-4,e-5]
%  Spacial phase profile modulation         Rminimum=[0.05]
%  Spacial frequencise                            vx=[6-:]
%  Slope and divergence of the beam            alpha=[1-5]
%
%  Offset=[]
%  BeamDecentering=[]
% Federico Vismarra 26/10/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=828;%CCD parameters
N=1485;
dx=0.0044; %mm one pixel is 4.4um
dy=0.0044; %mm one pixel is 4.4um
dz=0.1;
x=-N/2*dx:dx:N/2*dx-dx;
y=-M/2*dy:dy:M/2*dy-dy;
z=-1:dz:1; %arbitrary propagation length in 2mm?
[Xp,Yp]=meshgrid(x,y);
[Y,X,Z]=meshgrid(y,x,z);  %tensorial space for n(x,y,z)


%refractive index profile creation:
R0=0.1;  %minimum radius(mm)
alpha=2;  %slope
IT=exp(-X/alpha)+R0; 
Rinv=IT.^-1;
Dn0=0.0001;
Dn_xyz=Dn0*sqrt(pi*Rinv.^2).*exp(-(Z.^2+Y.^2).*Rinv.^2).^15.*(0.2+rand(size(X))); % più è sottile la differenza con l'aria meglio si vedono le modulazioni 


%la fase è una somma di n lungo z per ogni retta x,y
for i=1:length(x)
    for j=1:length(y)
phase(j,i)=sum(Dn_xyz(i,j,:))*2*pi/(685*10^-6)*dz; %everything set in mm 
%inverto gli indici perchè nel tensore le cose sono al contrario!
    end
end
phase(300:500,600:800)=phase(300:500,600:800)/2;
figure(1)
mesh(x,y,phase);
xlabel('X(mm)')
ylabel('Y(mm)')
zlabel('Phase(a.u.)')
title("Variation of the Phase_ Fronte d'onda Simulato");


%% Definisco il campo incidente tot e produco interferenza fringes
E0=10;
vx=7;

offset=2*pi-0.1;
Itot=(E0*cos(2*pi*vx*Xp+phase+offset)+2*E0).*exp(-(Xp.^2+Yp.^2)/(2*0.8^2));

figure(3);
mesh(Xp,Yp,Itot);
xlabel('Xp(mm)');
ylabel('Yp(mm)');
zlabel('I(arb)');

%% ANALISI ROW BY ROW OF ITOT
Phase=zeros(M,N); 
for i=1:M
R=PhaseRecovery2(Itot(i,:),x,dx,N);  % You can optimize this object playing with inner parameters
Phase(i,:)=R;
end
figure(3);
PHASE=Phase; %mind the offset have to be defined  and manually adjust... ( We can think a way to optimize it automatically)
% The constant slope is a weird effect of the quantization ???


% restore the flatness____AUTOMATIC SLOPE REMOVAL. It sample the slope at the edges!
select=1.5;
m=(PHASE(round(M/2-select/dy),round(N/2+select/dx))-PHASE(round(M/2-select/dy),round(N/2-select/dx)))/(x(round(N/2+select/dx))-x(round(N/2-select/dx)))
PHASE=PHASE-m*Xp;    



h=(PHASE(round(M/2-select/dy),round(N/2+select/dx))+PHASE(round(M/2-select/dy),round(N/2-select/dx)))/2;
PHASE=PHASE-h;




PHASE(:,1:round(N/2-2.5/dx))=0;

PHASE(:,round(N/2+2.5/dx):end)=0;

mesh(Xp,Yp,PHASE); % You need a slope adjustament WHY the constant term is the distance between the center of x=0 and the maximum?
xlabel('Xp(mm)');
ylabel('Yp(mm)');
zlabel('Phase Reconstructed(arb)');


