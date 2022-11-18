%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       **** Strong Field Electron XUV Interaction Simulator ****     %%%
%%%                                                                     %%%
%%%                                 MAIN 1D                             %%%
%%%                                                                     %%%
%%%                      Federico Vismarra, 12 May 2020                 %%%
%%% IF USED CITE: F Vismarra et al 2022 J. Phys. Photonics 4 034006     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               (Politecnico of Milano, Attosecond Research Center)             
% Special thanks to  M.Lucchini, Y.Wu and B.Moio
%
%  LEGEND: 
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           THINGS YOU CAN CUSTOMIZE   %
%                                      %
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  COMMENT:
%  SFEXIS is a simulation software which implements Strong Field Integral
%  to evaluate the electron ionization probability toward a free state.
%  The core of the process is a multi-quantum path interference (Propagators)
%  
%  The main assumpitions, that the software make, are STRONG FIELD APPROXIMATION (SFA) and
%  Single Electron Approximation (SEA) which neglects many-electron
%  correlation.
%  The software assume at the moment a flat fields profile for both XUV and IR (1D), this
%  approximation will be soon removed including 3D profile and angular
%  electron distribution.
%
% FEATURES OF SFEXIS:
% - GENERATION OF AN ABITRARY Gaussian XUV IR FIELD ( GVD TOD control)
% - RABBIT SIMULATION (ATP)
% - STREAKING SIMULATION (SAP)
% - SIDEBAND SIMULATION (SAP)
% - DYNAMIC DIPOLE EFFECT   
% - INTERACTION Target Selectivity (Ar,Ne,He)
% - LINEARLY POLARIZED FIELDS
% INFORMATION:
% Atomic units are used in the code
% While INPUTS are in more familiar units fs, nm, eV.
%
%
% 
% 
%
% SFEXIS by F.Vismarra
%  
%
%
clear all;
clc;
close all;
disp("********************************************");
disp("*******Welcome in SFEXISS v.ALPHA 1.0*******");
disp("**********************************F.V.******");
disp("**************MAIN 1D***********************")
disp("...Press anything to start...")
pause();
%% **************** CONSOLE *************** %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *** CONVERSION TO ATOMIC UNITS *** 
timeau=2.418884326505e-2; %(fs to a.u.)
c0= 2.99792458e8;         %(m/s)
EAu= 27.21138624;         %(au to eV)

disp('Starting Initialization...');          
       disp(" ");                                   
       disp("OPTIONS:");   
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.01; %Time Resolution (fs) (5 as is optimal ) (10 as for faster calculation)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
% Auxiliary Time
t=-5500+dt/2:dt:5500-dt/2;     %Time axis (fs)  ALWAYS EVEN

t_au=t/(timeau); % Time Axis (au)
dt_au=dt/timeau; % Resolution (au)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                SETTINGS                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IR Field                                                                                                %
CustomIR=0;              % (0) the IR is generated as a simple Gaussian pulse, (1) must be given in input.%
%Case 1
fileIR='';               %                                                                                %
%Case 0
FWHM_IR=20;              % IR Duration  (fs)                                                              %
lambda_IR=800;           % Wavelength IR in (nm)                                                          %
GDD_IR=0;               % Second Order of Dispersion  (fs^2)                                             %
TOD_IR=0;                % Third Order of Dispersion   (fs^3)                                             %
IR_CEP=0;                % (0-2*pi) shift with respect the center                                         %
Ipeak=1e11;              % Peak Intensity in W/cm^2                                                       %                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Selected lambaIR=%d (nm) FWHM_IR=%d (fs) GDD_IR=%d (fs^2) TOD_IR=%d (fs^2)',lambda_IR,FWHM_IR,GDD_IR,TOD_IR)
disp(" ")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XUV field (SAP) single pulse or ATP                                                                       %
% IR Field                                                                                                  %
CustomXUV=0;              % (0) the IR is generated as a simple Gaussian pulse, (1) must be given in input. %
%Case 1                                                                                                     %
fileXUV='';    %                                                                                              %
%Case 0%                                                                                                    %
FWHM_XUV=0.1;             % Durata XUV (fs)                                                                 %
xuvHH=21;                 % Central Harmonics for HHG                                                       %
GDD_XUV=0;                % Second Order of Dispersion  (fs^2)                                              %
TOD_XUV=0;                % Third Order of Dispersion   (fs^3)                                              %
XUV_CEP=0;                % (0-2*pi) shift with respect the center                                          %
% The Train Option Enable Train Of Pulses for Rabbit                                                        %  
TrainFlag=1;              % Rabbit experiment =1                                                            %
Extrapulses=5;            % x2 Number of pulses in the ATP                                                  %
odd=1;                    % Change the shape of the train                                                   %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Selected XUVHHG=%d FWHM_XUV=%d (fs) GDD_XUV=%d (fs^2) TOD_XUV=%d (fs^2)',xuvHH,FWHM_XUV,GDD_XUV,TOD_XUV)
disp(" ")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delay tau ( Move IR with respect XUV)                           %
FACTOR=10;                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dtau=FACTOR*dt;           % Resolution in delay axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            **** Atomic system definition ****                                 %
Atomic=1;               %(1) Ar: 15.75962 %(2) Ne: 21.5646 %(3) He: 24.58 (eV)  %                                                                                %
EnergeticRange=[5,30]; % in (eV)                                                %
ResE=0.05;              % in (eV)                                                %
flagdipole=0;          % 0 no dipole; 1 dipole                                  %
custom=0;              % 0 simple dipole WF; 1 custom dipole                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_fin=EnergeticRange(1):ResE:EnergeticRange(2);
v_au=sqrt(2*E_fin/EAu);
XuvEn=1.2398/(lambda_IR*10^(-3))*xuvHH;
fprintf("You are exiting with an XUV field with energy %d eV\n",XuvEn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************* SAVING OPTIONS ********************  %
CurrentDirectory=cd;                                 % 
Date=datetime                                        %
Date=datestr(Date);                                  % 
Date(Date==' ')=[];                                  %
Date(Date==':')='_';                                 %
SavingDirectory=[cd,'\Simulations'];                %
mkdir([SavingDirectory,'\' Date,'TEST1Ar']);               %
SavingDirectory=[cd,'\Simulations\', Date,'TEST1Ar'];     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(" ");
disp('XUV IR Fields Generation:');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLORBAR RESOLUTION
m=4; %Setting the resolution (m=1 Defeault)             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('jet_white.mat');
load('ePie_map.mat');
m=64*m;
cindex = linspace(1,64,m);
r = interp1([1:64],jet_white(:,1),cindex);
g = interp1([1:64],jet_white(:,2),cindex);
b = interp1([1:64],jet_white(:,3),cindex);
jet_white= [r' g' b'];


r = interp1([1:64],ePie_map(:,1),cindex);
g = interp1([1:64],ePie_map(:,2),cindex);
b = interp1([1:64],ePie_map(:,3),cindex);
ePie_map= [r' g' b'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% **************** INITIALIZATION *************** %%

%*******************IR*******************************% 
if CustomIR==0
        freq_IR=c0/(lambda_IR*10^(-9))*10^(-15);       % 1/fs
        freqIR_au=freq_IR*timeau;                      % Frequency IR in (1/a.u.)
        FWHMau_IR=FWHM_IR/(timeau);                    % Temporal width IR (electric field)  
        GDDau_IR=GDD_IR/(timeau^2);                    % Second Order of Dispersion (au^2)
        TODau_IR=TOD_IR/(timeau^3);                    % Third Order of Dispersion  (au^3)
        E0=sqrt(Ipeak/(3.51*10^16));                   % Conversion in atomic units of E0   I=0.5*E^2 (a.u.)
        % IR field generation 
        [A_IR,E_IR]=EField1D(FWHMau_IR,freqIR_au,GDDau_IR,TODau_IR,IR_CEP,t_au,dt_au);  % Normalized unit
                    
        % Conversion in atomic units

elseif CustomIR==1
    %Lo puoi caricare lungo quanto vuoi ATTENZIONE A CALCOLARE A_IR
    %CORRETTAMENTE DA E_IR 
    %E_IR E0=sqrt(Ipeak/(3.51*10^16)); % Conversion in atomic units of E0   
    %A_IR
%     % A resize
%         Asim_IR=zeros(length(tSim),1);
%         Asim_IR=interp1(t,A_IR,tSim);
% 
%         % E resize
%         Esim_IR=zeros(length(tSim),1);
%         Esim_IR=interp1(t,E_IR,tSim);
% 
%         % Check error on A area
%         error=trapz(Asim_IR);


end


  
% Resize of the vector ( A longer vector was defined for zero convergence purpouses)
flag=0;
nSTOP=2;

while nSTOP<7 && flag==0
% Check on the area confidence

tSim=-nSTOP*FWHM_IR+dt/2:dt:nSTOP*FWHM_IR-dt/2; % Effective Simulation Time (fs)!

% A resize
Asim_IR=zeros(length(tSim),1);
Asim_IR=E0*interp1(t,A_IR,tSim); % Normalization to the right field amplitude (a.u.)

% E resize
Esim_IR=zeros(length(tSim),1);
Esim_IR=E0*interp1(t,E_IR,tSim); % Normalization to the right field amplitude (a.u.)

% Check error on A area
error=trapz(Asim_IR);
if error<1e-6
    flag=1;
end 
nSTOP=nSTOP+1;


end
fprintf('A_IR area: %d\n',error)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*********************************XUV************** %
if CustomXUV==0;                                     
    freq_XUV=xuvHH*c0/(lambda_IR*10^(-9))*10^(-15);  % 1/fs
    freqXUV_au=freq_XUV*timeau;                      % Frequency IR in (1/a.u.)
    FWHMau_XUV=FWHM_XUV/(timeau);                    % Temporal width IR (electric field)  
    GDDau_XUV=GDD_XUV/(timeau^2);                    % Second Order of Dispersion (au^2)
    TODau_XUV=TOD_XUV/(timeau^3);                    % Third Order of Dispersion  (au^3)

    if TrainFlag==0
           [A_XUV,E_XUV]=EField1D(FWHMau_XUV,freqXUV_au,GDDau_XUV,TODau_XUV,XUV_CEP,t_au,dt_au); 

        elseif TrainFlag==1
            disp("RABBIT EXPERIMENT SELECTED");
           [A_XUV,E_XUV]=EField1D_train(FWHMau_XUV,freqXUV_au,GDDau_XUV,TODau_XUV,XUV_CEP,t_au,dt_au,freqIR_au,Extrapulses,odd,E_IR); 
    if odd==1 % Not Optimal Choice
        [E_IRupper,E_IRlower] = envelope(E_IR);
         E_XUV=E_XUV.*E_IRupper'/(max(E_IRupper)); 
       clear E_IRupper E_IRlower
       % THIS APPROXIMATION SHOULD BE IMPROVED SEE odd==0 om Efield1D_train
    end



        else
            disp("Error TrainFlag must be 0 or 1");
            return;

    end
elseif CustomXUV==1
    %dd
    
    
end


    % Rescaling and XUV reshaping
    Esim_XUV=interp1(t,E_XUV,tSim);
    Asim_XUV=interp1(t,A_XUV,tSim);
    A_area=trapz(Esim_XUV);
    fprintf('E_XUV area: %d\n',A_area)
    disp(" ");

disp('Field Genereted...I m plotting the fields');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Plot of the fields:
figure(1)
set(gcf,'Position',[50 50 800 700]);
subplot(2,1,2)
plot(tSim,Esim_XUV,'LineWidth',1,'Color',[0 0 1]);
title(['EXUV Pulse',' FWHM=', num2str(FWHM_XUV),' (fs)']);
xlabel('Time (fs)')
ylabel('Normalized Amplitude');
axis([min(tSim), max(tSim), -1, 1])
set(gca,'fontsize', 18)
subplot(2,1,1)
plot(tSim,Asim_IR,'LineWidth',3,'Color',[1 0 0]);
title(['AIR Pulse \lambda=' num2str(lambda_IR) ' (nm)' 'FWHM=' num2str(FWHM_IR) ' (fs)']);
xlabel('Time (fs)')
ylabel('Amplitude (a.u.)');
axis([min(tSim), max(tSim), -max(Asim_IR), max(Asim_IR)])
set(gca,'fontsize', 18)
cd(SavingDirectory);
saveas(figure(1),'XUVIRField.fig');
cd(CurrentDirectory);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*********************************DIPOLE ELEMENT**************************************************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(" ");
disp('...Quantum system initialiation...');
if flagdipole==1
    disp("OPTION: Dynamic Dipole Evaluation Selected (Hydrogenoid Model)");
end
switch Atomic
    
    case 1
        disp('You selected Ar(3p)');
        Ip=15.75962; %(eV)
        N=3; M=0; L=1 ; %Assumption on the interacting system
        Zeff=7;
    case 2
        disp('You selected Ne(2p)... You richy one!');
        Ip=21.5646; %(eV)
        N=2; M=0; L=1;
        Zeff=6;
    case 3
        disp('You selected He(1s).');
        Ip=24.5873; %(eV)
        N=1; M=0; L=0;
        Zeff=2;
    otherwise
        disp('Error! Unknown...I m setting He for you'); 
        Atomic=3;
        Ip=24.5873; %(eV)
        N=1; M=0; L=0;
        Zeff=2;
end


% STATIC DIPOLE
if flagdipole==1 && custom==0    
        % WAVEFUNCTION CALCULATION
        a=1/Zeff;  % Bohr radius (a.u.)

        % Angular part (Condon-Shortley)
        SphericalYlm = @(l, m, theta, phi) (-1)^m * sqrt((2 * l + 1) / (4 * pi) * ...
            factorial(l - abs(m)) / factorial(l + abs(m))) * ...
            AssociatedLegendre(l, m, cos(theta)) .* exp(1i * m * phi);

        Y = @(l, m, theta, phi) SphericalYlm(l, m, theta, phi);

        % Radial part
        R = @(n, l, r) sqrt((2 / (a * n))^3 * factorial(n - l - 1) / (2 * n * factorial(n + l))) .* ...
            exp(-r / (a * n)) .* (2 * r / (a * n)).^l * 1 / factorial(n - l - 1 + 2 * l + 1) .* ...
            AssociatedLaguerre(n - l - 1, 2 * l + 1, 2 * r / (a * n));

        % wave function
        psi = @(n, l, m, r, theta, phi) R(n, l, r) .* Y(l, m, theta, phi);

        %  setting the grid (OPTIMIZED)
        maxLength= 20;  % a.u.
        accuracy = 200;
        dim = linspace(-maxLength, maxLength, accuracy);
        [x, y, z] = meshgrid(dim, dim, dim);

        r = sqrt(x.^2 + y.^2 + z.^2);
        theta = acos(z ./ r);
        phi = atan2(y, x);
        PSI=psi(N, L, 0, r, theta, phi);
        thresholdOuter=0.01;
        thresholdInner=0.2;
        rBounding=3;


               F=figure;
               plotPsi(PSI,dim,thresholdOuter, thresholdInner, rBounding);
               title('Hydrogenoid state ');
        xlabel("x-axis (a.u.)");
        ylabel("y-axis (a.u.)");
        zlabel("z-axis (a.u.)");
              cd(SavingDirectory);
              saveas(F,'WaveFunction.fig');
              cd(CurrentDirectory);

        % STATIC DIPOLE CALCULATION (1D projection LINEARLY POLARIZED FIELDS)
        RES=10; %Parameter for performance increase
        FFTOPT=4;
        paux=linspace(-v_au(UR)*FFTOPT,v_au(UR)*FFTOPT,UR*RES); %Create an evenly spaced auxiliary vector;
        % Some manipulation for fft
        U=length(paux);
        dpaux=paux(2)-paux(1);
        paux=paux(1)-dpaux/2:dpaux:paux(U)-dpaux/2;
        U=length(paux);
        dipolestat1D=DipolestatCalculation1D(paux,PSI,x,y,z,dim);   



elseif flagdipole==1 && custom==1 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load("Ar_DipoleT_Khe.mat")                                              %
        hold on                                                                 %
        v_auimported=sqrt(2*EnergyAr/EAu);                                      %
        UR=length(v_auimported);                                                %
                                                                                %
        EnergyExtended=linspace(EnergyAr(1),EnergyAr(end),1e3)';                %
        v_auextended=sqrt(2*EnergyExtended/EAu);                                %
        Phase_Ar=smooth(interp1(EnergyAr,Phase_Ar,EnergyExtended,'pchip'),25);  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dipolestat=1.*exp(1i*Phase_Ar/(2*pi));   % ALWAYS CHECK% 
      
        figure(10)
        title('TimeDelay')
        plot(EnergyExtended(2:end),diff(Phase_Ar/(2*pi))/((EnergyExtended(2)-EnergyExtended(1))/EAu)*25)    
        xlabel('Energy (eV)');
        ylabel('Delay (as)');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Some manipulation for fft
        dipolestat1D=[flipud(dipolestat);dipolestat];
        paux=[-flipud( v_auextended); v_auextended];
        U=length(paux);

        
        
        
elseif flagdipole==0
        dipolestat1D=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Computation
%%%%%%%%%%%%%%%%%%%%%
% FIELD GENERATION  %
%%%%%%%%%%%%%%%%%%%%%


disp(" ");
disp('Creating the delayed signal matrix');
tminD=tSim(min(find(tSim>-(nSTOP-2)*FWHM_IR)));
tmaxD=tSim(min(find(tSim>(nSTOP-2)*FWHM_IR)));


delay=tminD:dtau:tmaxD; %(fs)

E_XUV_trace=zeros(length(delay),length(Esim_XUV)); % XUV IS THE ONE THAT IS MOVED!!!!
A_IR_trace=ones(length(delay),length(Asim_IR)).*Asim_IR;

for k=1:length(delay)
shiftAmount=(-length(delay)/2+k)*(dtau/dt);


E_XUV_trace(k,:)=E_XUV_trace(k,:)+circshift(Esim_XUV,shiftAmount);
if shiftAmount<0
E_XUV_trace(k,end+shiftAmount:end)=0;

elseif shiftAmount>0
E_XUV_trace(k,1:shiftAmount)=0;    
        
    else
        0;
end % Comment if you add Optional
end
%OPTIONAL FOR ENTERTEINMENT PU
%  subplot(2,1,2)
%  plot(tSim,E_XUV_trace(k,:));
%  title('XUV Pulse');
%  xlabel('Time (fs)')
%  ylabel('Normalized Amplitude');
%  axis([min(tSim), max(tSim), -1, 1])  
% 
%  subplot(2,1,1)
%  plot(tSim,A_IR_trace(k,:))
%  title('IR Pulse');
%  xlabel('Time (fs)')
%  ylabel('Amplitude (a.u.)');
%  axis([min(tSim), max(tSim), -max(Asim_IR), max(Asim_IR)])
%  
%  pause(0.001);
%  end
 clear A_IR A_XUV E_IR E_XUV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************** FINAL INTEGRATION **********************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1;
wantedTime=max(tSim);                             % Final at time at which you look the effect of the interaction
tmin=min(tSim);
tmax=max(tSim);

% INITIALIZATION
av=zeros(length(delay),length(v_au));                 % Velocity Transition amplitude in time.
WE=zeros(length(delay),length(v_au));                 % Energy   Transition amplitude in time.

 time=tmin:dt:wantedTime;
 EatomicReal=zeros(length(delay),length(time),length(v_au));
 flag=1;
%      dip=ones(length(delay),length(time),length(v_au))*1i; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******      % DIPOLE FEATURE CALCULATION %   *********%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flagdipole==1       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        dipoleAux=ones(length(time),U)*1i;
        dipoleA=ones(length(time),length(v_au))*1i;
        for j=1:length(time)
        dipoleAux(j,:)=dipoleShift(dipolestat1D,Asim_IR(j),paux); 
        dipoleA(j,:) = interp1(paux,dipoleAux(j,:),v_au);
%         figure(7)
%         plot(v_au,unwrap(angle(dipoleA(j,:))),'ob')
%         hold on
%         plot(v_au,abs(dipoleA(j,:)),'or');
%         pause(0.001)
%         hold off
        end
        disp('Dynamic Dipole Calculated........')
        disp(' ');
           
%         for b=1:length(delay)
%         dip(b,:,:)=dipoleA; %used as a function of delay and time.
%         b
%         end
          disp("Dipole Calculated");
          
                                    fp=abs(dipoleA); % PLOT DYNAMIC DIPOLE
                                    fig_trace=figure(2);
                                    set(gcf,'Position',[100 100 800 800])
                                    pcolor(time,E_fin,fp'); % x y z
                                    ylabel('Energy state (eV)');
                                    xlabel('Time (fs)');
                                    set(gca, 'Layer', 'top');
                                    shading flat;
                                    colormap(ePie_map);
                                    cb=colorbar;
cd(SavingDirectory);
saveas(fig_trace,['DynamicDipole','.fig']);
cd(CurrentDirectory);
                   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START VELOCITY SCAN: REPEAT FOR ALL THE VELOCITIES%%

Saux=zeros(length(delay),length(time));                  % Semiclassical Action.
EXUVfin=ones(length(delay),length(time));                % Final Field
T=ones(length(delay),length(time)).*(time/(timeau));     % Time

EXUVfin=EXUVfin.*E_XUV_trace(:,1:length(time));          % XUV TRACE COMPOSITION


for m=1:length(v_au)
    tic;
    % You can customize this cycle if you want to speed up the
    % computation.. large memory needed, though.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******* Semiclassical action calculation ******%
arg=zeros(length(delay),length(time)); 

tempIR=A_IR_trace(:,1:length(time));
tempIR2=(-tempIR(:,k)+tempIR);

Saux=tempIR2+v_au(m);
Saux=0.5*flip(cumsum(flip(Saux.^2,2),2),2)*dt/(timeau);
clear tempIR tempIR2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Final CALCULATION %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flagdipole==1
dip=ones(length(delay),length(time)).*dipoleA(:,m)';

else
    dip=1;
end


arg=dip.*exp(-1i*Saux).*exp(1i*Ip/(EAu)*(T-tmin/(timeau))).*EXUVfin+arg;
av(:,m)=-1i*sum(arg,2)*(dt/timeau); % Alternatively you can use trapz("")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MEMORY INFORMATIONcxzzxc
Info=whos;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Saux arg;

C=toc;
disp(" ")

fprintf("Time Evaluation VSCAN n=%d/%d time take",m,length(v_au),C)
disp(" ");
WE(:,m)=abs(av(:,m)).^2*v_au(m); % av to aE: Conversion to real energy MAP
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING

disp("I'm Plotting and Saving the FINAL Result...");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    WE TOT IR  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOT WE TOT WITH IR PROFILE
fp=WE;
fp=fp/max(max(fp));
fig_trace=figure(3);
set(gcf,'Position',[100 100 800 800])
pcolor(delay,E_fin,fp'); % x y z
ylabel('Electron Energy (eV)');
xlabel('Delay (fs)');
set(gca, 'Layer', 'top');
shading flat;
caxis([0 1]);
colormap(jet_white);
cb=colorbar;

axis([tminD,tmaxD,EnergeticRange(1),EnergeticRange(2)]);
title(['XUV Duration= ' num2str(FWHM_XUV) '(fs) IR Duration=' num2str(FWHM_IR) '(fs) Ip=' num2str(Ip),'(eV) HH=' num2str(xuvHH), 'Ipeak=' num2str(Ipeak/10^12) '(TW/cm2)', ' FINAL PLOT'  ])
cd(SavingDirectory);
saveas(fig_trace,['WeFINALE','.fig']);
cd(CurrentDirectory);



%SAVING DATA:
cd(SavingDirectory);
save(['Wetot.mat'],'WE','fp','delay','E_fin')
cd(CurrentDirectory);
