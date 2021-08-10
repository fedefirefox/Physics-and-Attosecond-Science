clear all
clc
close all
%% Parameters definition
timeau=2.418884326505e-2; %(fs to a.u.)
c0= 2.99792458e8*10^(-15)*10^(9);         %(nm/fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.01; %Time Resolution (fs) (5 as is optimal ) (10 as for faster calculation)%
lambda_IR=800;           %Central frequency  (nm)                                    %
FWHM_IR=10;              % IR Duration  (fs)                                    %
GDD_IR=10;                % Second Order of Dispersion  (fs^2)                   %
TOD_IR=0;                % Third Order of Dispersion   (fs^3)                   %
IR_CEP=0;                % (0-2*pi) shift with respect the center               %
Ex_0=1;                    % Initial field 
Ey_0=0;                    % Initial field
phase=0;                   % Initial phase between the x-component and y-component (polarization)
thetaA=45;                %Angle between the x-component and the ordinary axis of the thick plate
theta=thetaA/(360)*2*pi;
gammaA=0;                %Angle between the x-component and the ordinary axis of the QWP
gamma=gammaA/(360)*2*pi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_c=2*pi/(lambda_IR);

% Auxiliary Time
t=-5500:dt:5500;     %Time axis (fs)  ALWAYS EVEN
t_au=t/(timeau); % Time Axis (au)
dt_au=dt/timeau; % Resolution (au)
f_au=1/timeau;


%% Electric field calculation temporal profile

        freq_IR=c0/(lambda_IR);       % 1/fs
        freqIR_au=freq_IR*timeau;                      % Frequency IR in (1/a.u.)
        FWHMau_IR=FWHM_IR/(timeau);                    % Temporal width IR (electric field)  
        GDDau_IR=GDD_IR/(timeau^2);                    % Second Order of Dispersion (au^2)
        TODau_IR=TOD_IR/(timeau^3);                    % Third Order of Dispersion  (au^3)
        
        % IR field generation 
        [E_IR,A_f,fau]=EField1D(FWHMau_IR,freqIR_au,GDDau_IR,TODau_IR,IR_CEP,t_au,dt_au);  % Normalized unit        
        f=fau*f_au; %Frequency axis in PHz
        figure(1)
        subplot(2,1,1)
        plot(t,E_IR) %this is the temporal shape common to both pol. at the beginning.
        axis([-30 30 -1 1])
        subplot(2,1,2)
        yyaxis left
        plot(f,abs(A_f))
        axis([0 1 0 1])
        yyaxis right
        plot(f,-unwrap(angle(A_f)))
 %%
 %Ideal Retarder

Eo_0=Ex_0*exp(1i*phase)*cos(theta)*A_f+Ey_0*sin(theta)*A_f;  % Initial condition of field among the ordinary axis
Ee_0=-Ex_0*exp(1i*phase)*sin(theta)*A_f+Ey_0*cos(theta)*A_f; % Initial condition of field among the extraordinary axis

Ex_f0=Ex_0*exp(1i*phase)*A_f;
Ey_f0=Ey_0*A_f;



n0_f=ones(length(f),1)'; %ordinary axis for the given material (slow-axis) e.g. glass
deltaT=9; %Time delay between the two components
deltaL=c0*deltaT;
PHI=2*pi*deltaT*(f-freq_IR);%(ones(length(f),1)*pi)'; %Data taken from BhalleWaveplate 1) HWP 2) QWP 3) Retarder
z0=2*10^4; %Thickness of the object

z=2*10^4; %100 um thick object
delta=PHI/z0;


D=exp(-1i*(2*pi*f.*n0_f/c0-k_c)*z); %Dispersion term common to both
Eo=Eo_0*exp(-1i*k_c*z).*D;  
Ee=Ee_0*exp(-1i*k_c*z).*D.*exp(1i*delta*z);



Ex=cos(theta)*Eo-Ee*sin(theta);
Ey=sin(theta)*Eo+cos(theta)*Ee;


Ex_t=fftshift(ifft((Ex)));
Ey_t=fftshift(ifft((Ey)));

Ex_0t=fftshift(ifft((Ex_f0)));
Ey_0t=fftshift(ifft((Ey_f0)));

normF=max(sqrt(Ex_0t.^2+Ey_0t.^2));



Ex_0t=Ex_0t/normF;
Ey_0t=Ey_0t/normF;
Ex_t=Ex_t/normF;
Ey_t=Ey_t/normF;


figure(3)
plot(t,Ex_t,'r');
hold on
plot(t,Ey_t,'b');
hold off
axis([30 80 -1 1])


figure(4)

g1=plot3(t*c0,real(Ex_t),real(Ey_t))
hold on
g2=plot3(t*c0,real(Ex_0t),real(Ey_0t))
title('Propagation in a retarder')
xlabel ("space (nm)");
ylabel ("E_y ( arb.)");
zlabel ("E_x ( arb.)");
grid on
set(gca,'YTick',[-1 1], 'ZTick',[-1 1]);
set(g1,'LineWidth',1.5);
set(g2,'LineWidth',1.5);
axis([-30*c0 80*c0 -1 1 -1 1]) % caso 3D
hold off


%% QWP
Ex_f2=(fft(ifftshift(Ex_t)));
Ey_f2=(fft(ifftshift(Ey_t)));

Eo_0=Ex_f2*cos(gamma)+Ey_f2*sin(gamma);  % Initial condition of field among the ordinary axis
Ee_0=-Ex_f2*sin(gamma)+Ey_f2*cos(gamma); % Initial condition of field among the extraordinary axis

n0_f=ones(length(f),1)'; %ordinary axis for the given material (slow-axis) e.g. glass

PHI=(ones(length(f),1)*pi/2)'; %Data taken from BhalleWaveplate 1) HWP 2) QWP 3) Retarder
z0=2*10^4; %Thickness of the object

z=2*10^4; %100 um thick object
delta=PHI/z0;


D=exp(-1i*(2*pi*f.*n0_f/c0-k_c)*z); %Dispersion term common to both
Eo=Eo_0*exp(-1i*k_c*z).*D;  
Ee=Ee_0*exp(-1i*k_c*z).*D.*exp(1i*delta*z);



Ex=cos(theta)*Eo-Ee*sin(theta);
Ey=sin(theta)*Eo+cos(theta)*Ee;


Ex_t2=real(fftshift(ifft((Ex))));
Ey_t2=real(fftshift(ifft((Ey))));

Ex_0t=real(fftshift(ifft((Ex_f2))));
Ey_0t=real(fftshift(ifft((Ey_f2))));



normF=max(sqrt(Ex_0t.^2+Ey_0t.^2));



Ex_0t=Ex_0t/normF;
Ey_0t=Ey_0t/normF;
Ex_t2=Ex_t2/normF;
Ey_t2=Ey_t2/normF;


figure(7)
plot(t,Ex_t2,'r');
hold on
plot(t,Ey_t2,'b');
hold off
axis([110,155,-1,1])


figure(8)

g1=plot3(t*c0,real(Ex_t2),real(Ey_t2))
axis([110*c0 155*c0 -1 1 -1 1]) % caso 3D
title('Propagation in QWP')
xlabel ("space (nm)");
ylabel ("E_y ( arb.)");
zlabel ("E_x ( arb.)");
grid on
set(gca,'YTick',[-1 1], 'ZTick',[-1 1]);
set(g1,'LineWidth',1.5);
set(g2,'LineWidth',1.5);

hold off