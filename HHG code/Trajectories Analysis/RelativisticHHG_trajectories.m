%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HHG classical trajectories with magnetic field
%
% Federico Vismarra July 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all
%% Constants
% *** CONVERSION TO ATOMIC UNITS *** 
timeau=2.418884326505e-2; %(fs to a.u.)
c0= 2.99792458e8;         %(m/s)
EAu= 27.21138624;         %(au to eV)
xau=5.291e-2;              % (nm to au)
vau=2.18e6;               %(au to m/s);   
cau=c0/vau;               % (speed of light in atomic units)
%% Time definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.01; %Time Resolution (fs) (5 as is optimal ) (10 as for faster calculation)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
% Auxiliary Time (for fft)
t=-5500+dt/2:dt:5500-dt/2;     %Time axis (fs)  ALWAYS EVEN
t_au=t/(timeau); % Time Axis (au)
dt_au=dt/timeau; % Resolution (au)
clear t;

%% Input of IR pulse (in Fourier space) B amd E

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IR Field                                                                                                %
CustomIR=0;              % (0) the IR is generated as a simple Gaussian pulse, (1) must be given in input.%
%Case 1
fileIR='';               %                                                                                %
%Case 0
FWHM_IR=200;               % IR Duration  (fs)                                                              %
lambda_IR=3500;           % Wavelength IR in (nm)                                                           %
GDD_IR=0;                % Second Order of Dispersion  (fs^2)                                             %
TOD_IR=0;                % Third Order of Dispersion   (fs^3)                                             %
IR_CEP=0;               % (0-2*pi) shift with respect the center                                         %
Ipeak=1e15;              % Peak Intensity in W/cm^2                                                       %                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Selected lambaIR=%d (nm) FWHM_IR=%d (fs) GDD_IR=%d (fs^2) TOD_IR=%d (fs^2)',lambda_IR,FWHM_IR,GDD_IR,TOD_IR)
disp(" ")

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %            **** Atomic system definition ****                                                   %
% Atomic=1;             %(1) Ar: 15.75962 %(2) Ne: 21.5646 %(3) He: 24.58 (eV) %(4) Kr:13.9999 (eV) %                                                                                %
% EnergeticRange=[5,50];% in (eV)                                                                   %
% ResE=0.05;            % in (eV)                                                                   %
% flagdipole=0;         % 0 no dipole; 1 dipole                                                     %
% custom=1;             % 0 simple dipole WF; 1 custom dipole                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ip=15.75962; % (eV) 
%*******************IR*******************************% 
if CustomIR==0
        freq_IR=c0/(lambda_IR*10^(-9))*10^(-15);       % 1/fs
        freqIR_au=freq_IR*timeau;                      % Frequency IR in (1/a.u.)
        FWHMau_IR=FWHM_IR/(timeau);                    % Temporal width IR (electric field)  
        GDDau_IR=GDD_IR/(timeau^2);                    % Second Order of Dispersion (au^2)
        TODau_IR=TOD_IR/(timeau^3);                    % Third Order of Dispersion  (au^3)
        E0=sqrt(Ipeak/(3.51*10^16));                   % Conversion in atomic units of E0   I=0.5*E^2 (a.u.)
        % IR field generation 
        [A_IR,E_IRx]=EField1D(FWHMau_IR,freqIR_au,GDDau_IR,TODau_IR,IR_CEP,t_au,dt_au);  % Normalized unit
         E_IRx=E_IRx*E0;  
        [T_IR,B_IRy]=EField1D(FWHMau_IR,freqIR_au,GDDau_IR,TODau_IR,IR_CEP,t_au,dt_au); 
         B_IR_y=B_IRy/cau*E0;
      
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

Up=E0/(2*pi*freqIR_au) %ponderomotive energy in a.u.
Ecutoff=Ip+Up*3.17*EAu; %cut off energy (estimated)
opt_cycle = lambda_IR*1e-9/c0*1e15; % (fs)

%% Resolution of equation of motion
% 
% taux=-FWHMau_IR:dt_au:2*FWHMau_IR; % Time used in the resolution of the propagation equation
% temission=-FWHMau_IR:dt_au:2*FWHMau_IR; % same as before buy it explicitly counts the moment of emission
taux=0:dt_au:(2*opt_cycle)/timeau; % Time used in the resolution of the propagation equation
temission=0:dt_au:(2*opt_cycle)/timeau; % same as before buy it explicitly counts the moment of emission
 
cmap=jet(length(temission)); % Colormap

E_IR_au=interp1(t_au,E_IRx,taux); % Interpolation of E_IR
B_IR_au=interp1(t_au,B_IR_y,taux);
v0=[0;0];
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
figure(1)
for i=1:length(temission)
By=B_IR_au;
Ex=E_IR_au;
    
    for k=1:length(taux);
        
        if k<i
       vx(i,k)=nan;
       vz(i,k)=nan;
       v0x(i,k)=nan;
       v0z(i,k)=nan;
       By(k)=0;
       Ex(k)=0;
        else
        

        end
    end
        [tnum2,v]=ode113(@(t,v)ODEvsparticle_2(t,v,By,Ex,taux),taux,v0,opts);% Velocity in a.u.
            v0x(i,:)=v(:,1)';
            v0z(i,:)=v(:,2)';
     
            [tnum2,vrel]=ode113(@(t,v)ODEvsparticle_rel(t,v,By,Ex,taux,v0x(i,:),v0z(i,:)),taux,v0,opts);
            vx(i,:)=vrel(:,1);
            vz(i,:)=vrel(:,2);
    for k=1:length(taux);
        if k<i
        vx(i,k)=nan;
        vz(i,k)=nan;
        else
        end
        end
        
    plot(taux*timeau,0.5*(vx(i,:).^2+vz(i,:).^2)*EAu,'o','Color',cmap(i,:))
    hold on
    i
end
xlabel('Time (fs)')
ylabel('Kinetic energy (eV)')
title('Kinetic energy of all electrons')

figure(2)
for i=1:length(temission)

    
    for k=1:length(taux);
        
        if k<i
       xx(i,k)=nan;     
       zz(i,k)=nan;
        else
       xx(i,k)= trapz(vx(i,i:k))*dt_au;
       zz(i,k)= trapz(vz(i,i:k))*dt_au;
        end

       
            
    end
    subplot(2,1,1)
    yyaxis left
    plot(taux*timeau,xx(i,:)*xau,'-','Color',cmap(i,:))
    hold on
    subplot(2,1,2)
    plot(taux*timeau,zz(i,:)*xau,'-','Color',cmap(i,:))
    hold on
end
   subplot(2,1,1)
   yyaxis right
   plot(taux*timeau,E_IR_au)
   ylabel('E-field (a.u.)')
   yyaxis left
   plot(taux*timeau,zeros(length(taux),1),'k-')
   hold off
      ylabel('Space (nm)')
      xlabel('Time (fs)')
 subplot(2,1,2)
 yyaxis right
   plot(taux*timeau,E_IR_au)
   ylabel('E-field (a.u.)')
   yyaxis left
   plot(taux*timeau,zeros(length(taux),1),'k-')
   hold off
      ylabel('Space (nm)')
      xlabel('Time (fs)')
 
title('Trajectory of all electrons')



   %% Recombination and trajectories analysis 
   % here only the trajectories that recombined are considered when 
   
   vtot=sqrt(vx.^2+vz.^2);
   flag=0; %useful 
   %taux_enhanced=-FWHMau_IR:dt_au/10:2*FWHMau_IR;
   taux_enhanced=0:dt_au/10:(2*opt_cycle)/timeau;
   N=1;
   
for i=1:length(temission)
    xx_enhanced(i,:)=interp1(taux,xx(i,:),taux_enhanced);
    zz_enhanced(i,:)=interp1(taux,zz(i,:),taux_enhanced);
    v_enhanced(i,:)=interp1(taux,vtot(i,:),taux_enhanced);
    
 for k=1:length(taux_enhanced)
 
     if abs(xx_enhanced(i,k))>1.1 && flag==0   %abs(sqrt(xx_enhanced(i,k).^2+zz_enhanced(i,k).^2))>10.1 && flag==0 NB interesting.
     flag=1;
     end
     
     if flag==1 && abs(xx_enhanced(i,k))<1 %flag==1 && abs(sqrt(xx_enhanced(i,k).^2+zz_enhanced(i,k).^2))<10
     flag=2; %the electron has recombined, set at zero the other trajectories
     INDEX=k;
     N=N+1;   % recollision counts
     end
     
% Extra trajectories ..... 
%  if flag==2 && abs(xx_enhanced(i,k))>1
%     flag=3;   % Electron exit from atomic interaction a second time
%  end
% 
%  if flag==3 && abs(xx_enhanced(i,k))<0.2
%     flag=4;   % Electron re-encounter the atom a second time
%     INDEX=k;
%     Nlong=Nlong+1;
%     Nshort=Nshort-1;
%     
%  end
% 
%  if flag==4 && abs(xx_enhanced(i,k))>1
%     flag=5;   % Electron exit from atomic interaction a third time
%  end
% 
%  if flag==5 && abs(xx_enhanced(i,k))<0.5
%     flag=6;   % Electron re-encounter the atom a third time
%     Nlong=Nlong-1;
%     Nshort=Nshort-1;
%     NextraL=NextraL+1;
%  end

    if flag==2
    xx_enhanced(i,k)=nan; 
    zz_enhanced(i,k)=nan;
    v_enhanced(i,k)=nan;
    end
     
 end

 if flag==1
    xx_enhanced(i,:)=nan;
    zz_enhanced(i,:)=nan;
    v_enhanced(i,:)=nan;
 end
 
 if flag==3 || flag==5
    xx_enhanced(i,INDEX:end)=nan;
    zz_enhanced(i,INDEX:end)=nan;
    v_enhanced(i,INDEX:end)=nan;
 end
 flag=0;
 
end


figure(4)
for i=1:length(temission)  
%         yyaxis left
    subplot(2,1,1)
    plot(taux_enhanced*timeau,xx_enhanced(i,:)*xau,'-','Color',cmap(i,:))
    hold on
    subplot(2,1,2)
    plot(taux_enhanced*timeau,zz_enhanced(i,:)*xau,'-','Color',cmap(i,:))
    hold on
end
subplot(2,1,1)
title('Recombining trajectories in x')
xlabel('Time (fs)')
ylabel('Space (nm)')
subplot(2,1,2)
title('Recombining trajectories in z')
xlabel('Time (fs)')
ylabel('Space (nm)')
% end
% yyaxis right
%    plot(temission*timeau,E_IR_au)
%    yyaxis left
%    plot(temission*timeau,zeros(length(temission),1),'k-')
figure(5)

    for i=1:length(temission)
    plot(taux_enhanced*timeau,0.5*v_enhanced(i,:).^2*EAu,'o','Color',cmap(i,:))
    hold on
    end
 
title('Kinetic energy of recombining trajectories')
xlabel('Time (fs)')
ylabel('Kinetic energy of the electron (eV)')



 %% Kinetic energy of the recombining electron
 
 v_enhanced(isnan(v_enhanced))=0;
 
 % For each recombining trajectory find emission time, recombination time
 % and photon energy.
 
 for i=1:length(temission)
    K=find(abs(v_enhanced(i,:))>0);
    if isempty(K)
    index=1;
    else 
    index=K(end); %recombination time is the last non zero entry.             
    index2=K(1); %v_x of non recombining trajectories (or not yet started) is zero.
    recTime(i)=taux_enhanced(index);
    emTime(i)=taux_enhanced(index2);
    end

    vxF(i,:)=zeros(length(taux_enhanced),1);
    vxF(i,:)=nan;
    vxF(i,index)=v_enhanced(i,index);
    velocityD(i)=v_enhanced(i,index);
 
 end

 figure(12)
for i=1:length(temission)  
 INDEX=round(abs(0.5*velocityD(i).^2*EAu+Ip)/max(abs(0.5*velocityD.^2*EAu+Ip))*length(cmap));
 if INDEX==0
    INDEX=1;
 end
    subplot(2,1,1)
    plot(taux_enhanced*timeau,xx_enhanced(i,:)*xau,'Color',cmap(INDEX,:),'Linewidth',2)
    hold on
    subplot(2,1,2)
    plot(taux_enhanced*timeau,zz_enhanced(i,:)*xau,'Color',cmap(INDEX,:),'Linewidth',2)
    hold on
end
 subplot(2,1,1)
 title('Recombining trajectories in x')
 xlabel('Time (fs)')
 ylabel('Space (nm)')
 subplot(2,1,2)
  title('Recombining trajectories in z')
 xlabel('Time (fs)')
 ylabel('Space (nm)')
% %% Do graphs
%   figure(6)
%     for i=1:length(temission)
%  
%     plot(taux_enhanced*timeau,vxF(i,:).^2*EAu+Ip,'o','Color',cmap(i,:))
%     hold on
%  
%     end
% ylabel('Energy of emitted photon (eV)')    
%   yyaxis right
%   plot(taux*timeau,E_IR_au)
% ylabel('E-field amplitude (a.u.)')   
% title('Emitted Photon Energy')
% xlabel('Time (fs)')
%   
%   
%   
% 
%   figure(7)
%   hist(velocityD.^2*EAu+Ip,200)
%   
%   
%   figure(8)
%   plot(emTime*timeau,recTime*timeau,'o')
%   xlabel('Emission time (fs)')
%   ylabel('Recombination time (fs)')
%   
  figure(9)
  for i=1:length(velocityD)-1
  plot(emTime(i)*timeau,0.5*velocityD(i).^2*EAu+Ip,'o','Color',cmap(i,:))
  hold on
  end
  xlabel('Emission time (fs)')
  ylabel('Photon energy (eV)')
  
  
figure(10)
  for i=1:length(velocityD)-1
 INDEX=round(abs(0.5*velocityD(i).^2*EAu+Ip)/max(abs(0.5*velocityD.^2*EAu+Ip))*length(cmap));
 if INDEX==0
    INDEX=1;
 end
  plot(recTime(i)*timeau,0.5*velocityD(i).^2*EAu+Ip,'o','Color',cmap(INDEX,:),'Linewidth',2)
  hold on
  end
  xlabel('Recombination time (fs)')
  ylabel('Photon energy (eV)') 
title('Emitted Photon Energy')
xlabel('Recombination time (fs)')




%% select trajectory with brush in recombination time 
%

recTimefs=recTime*timeau;
figure
yyaxis left
plot(recTimefs, 0.5*velocityD(1:end-1).^2*EAu+Ip,'o')
ylabel('Photon energy (eV)')
hold on 
yyaxis right
plot(temission*timeau,E_IR_au)
ylabel('Normalized Electric field ')
xlabel('Recombination time (fs)')


figure
ppt1=plot(recTimefs, 0.5*velocityD(1:end-1).^2*EAu+Ip,'o')
xlabel('Recombination time (fs)')
ylabel('Photon energy (eV)') 

brush on  
pause
title('select trajectory with Brush')
ind_rec=find(get(ppt1,'BrushData'));
brush_rec = logical(get(ppt1,'BrushData'));
xd_rec=get(ppt1,'XData');
yd_rec=get(ppt1,'YData'); 
brushedx_rec=xd_rec(brush_rec);
brushedy_rec=yd_rec(brush_rec);
figure
plot(brushedx_rec,brushedy_rec,'o')
xlabel('Recombination time (fs)')
ylabel('Photon energy (eV)')
title('Selected Trajectory')
brush off


%% fit the trajectory

recTime_oc=brushedx_rec/(opt_cycle);
yVelocity=brushedy_rec;
xTime=recTime_oc ;
f=figure 
plot(xTime,yVelocity)
hold on 
yyaxis right
plot(taux*timeau/opt_cycle , E_IR_au/max(E_IR_au), '-', 'Color', cmap(i,:))
plot(taux*timeau/opt_cycle, 0.*E_IR_au/max(E_IR_au), 'Linestyle','--')
hold off

ylabel('Photon energy (eV)') 
xlabel('Recombination time optical cycle') 
title('Emitted Photon Energy, trajectory will be shifted of 1/4 of opt cycle')


figure 
answer=[];
while ~strcmp(answer,'Yes')
xTime_flip = flip(xTime);
xTime_shift=xTime_flip-0.25;
yVelocity_flip=flip(yVelocity);
plot(xTime_shift,yVelocity_flip)
xlabel('time')
ylabel('Velocity')


title('select degree of polinomial fit in the Promt')
prompt = ['select degree of polinomial fit '];
degree = input(prompt);
%close(f)


%perform a polyfit of the trajectory
p_degree = polyfit(xTime_shift,yVelocity_flip, degree);
%convert the polyfit into syms 
syms x 
fit_velocity(x)=poly2sym(p_degree)
curveIp(x)= Ip+0*x
curveMax(x)= max(yVelocity_flip)+0*x
figure
fplot(fit_velocity)
hold on 
plot(xTime_shift,yVelocity_flip,'o','LineWidth',3)
legend('fit', 'electron trajectory')
xlabel('recombination time')
ylabel('energye (eV)')
xlim([xTime_shift(1), xTime_shift(end)])
ylim([min(yVelocity_flip),max(yVelocity_flip)])
quest = 'Is The Polynomial Fit good?';
answer = questdlg(quest,'Yes','No')
end
%% calculate tc and tp SHORT

omega_C=max(yVelocity_flip); 
omega_P=Ip; 
mean_Omega= ((omega_C+omega_P)/2)+0*x
tmean=vpasolve(fit_velocity==mean_Omega, x,[xTime_shift(1),xTime_shift(find(yVelocity_flip==max(yVelocity_flip)))])
figure
hold on
fplot(mean_Omega)
fplot(fit_velocity)
plot(tmean,(omega_C+omega_P)/2,'*','LineWidth',5)

xlim([xTime_shift(1), xTime_shift(end)])
ylim([min(yVelocity_flip),max(yVelocity_flip)])

% %calculate the derivatives
assume(x, 'real')
der1(x)= (diff(fit_velocity, x, 1));
der2(x) = (diff(fit_velocity, x, 2));

%calculate the tangent 
x_tang_point= eval(tmean);
tangent =  der1(x_tang_point)*(x-x_tang_point) + fit_velocity(x_tang_point); 

figure(17)
hold on
fplot(tangent,'LineWidth',2)
fplot(fit_velocity,'LineWidth',2)
plot(xTime_shift, yVelocity_flip,'ko')
plot(xTime_shift, Ip.*ones(size(yVelocity_flip)), 'Linestyle','--')
plot(xTime_shift,  max(yVelocity_flip).*ones(size(yVelocity_flip)), 'Linestyle','--')
tp=vpasolve(tangent==curveIp, x)
tc=vpasolve(tangent== curveMax, x)
plot(tc,omega_C,'*','LineWidth',4)
plot(tp,omega_P,'*','LineWidth',4)
legend('tangent','fit', 'data','Ip', 'maxCurve','tc','tp')
xlim([xTime_shift(1), xTime_shift(end)])
ylim([Ip-0.5,max(yVelocity_flip)+0.5])
xlabel('Recombination Time (optical cycle)')
ylabel('XUV photon energy (eV)')



%% calculate tc and tp LONG

omega_C=max(yVelocity_flip); 
omega_P=Ip; 
mean_Omega= ((omega_C+omega_P)/2)+0*x
tmean=vpasolve(fit_velocity==mean_Omega, x,[xTime_shift(find(yVelocity_flip==max(yVelocity_flip))),xTime_shift(end)])
figure
hold on
fplot(mean_Omega)
fplot(fit_velocity)
plot(tmean,(omega_C+omega_P)/2,'*','LineWidth',5)

xlim([xTime_shift(1), xTime_shift(end)])
ylim([min(yVelocity_flip),max(yVelocity_flip)])

% %calculate the derivatives
assume(x, 'real')
der1(x)= (diff(fit_velocity, x, 1));
der2(x) = (diff(fit_velocity, x, 2));

%calculate the tangent 
x_tang_point= eval(tmean);
tangent =  der1(x_tang_point)*(x-x_tang_point) + fit_velocity(x_tang_point); 

figure (18)
hold on
fplot(tangent,'LineWidth',2)
fplot(fit_velocity,'LineWidth',2)
plot(xTime_shift, yVelocity_flip,'ko')
plot(xTime_shift, Ip.*ones(size(yVelocity_flip)), 'Linestyle','--')
plot(xTime_shift,  max(yVelocity_flip).*ones(size(yVelocity_flip)), 'Linestyle','--')
tp=vpasolve(tangent==curveIp, x)
tc=vpasolve(tangent== curveMax, x)
plot(tc,omega_C,'*','LineWidth',4)
plot(tp,omega_P,'*','LineWidth',4)
legend('tangent','fit', 'data','Ip', 'maxCurve','tc','tp')
xlim([xTime_shift(1), xTime_shift(end)])
ylim([Ip-0.5,max(yVelocity_flip)+0.5])
xlabel('Recombination Time (optical cycle)')
ylabel('XUV photon energy (eV)')


%% cut off law
omega_C_au= omega_C/EAu 
a=3.17/(4*4*pi^2);

omega_C_ref= a*(E0*opt_cycle/timeau)^2 + Ip/EAu


%%
function dvdt = ODEvsparticle_2(t,v,By,Ex,tIN)
Byt = interp1(tIN,By,t); % Interpolate the data set (gt,g) at time t
Ext = interp1(tIN,Ex,t); % Interpolate the data set (gt,g) at time t

dvdt=[-Ext+v(2).*Byt;-v(1).*Byt]; % 
end


function dvdt = ODEvsparticle_rel(t,v,By,Ex,tIN,v0x,v0z)
vau=2.18e6;               %(au to m/s);   
c0= 2.99792458e8;         %(m/s)
cau=c0/vau;               % (speed of light in atomic units)
Byt = interp1(tIN,By,t); % Interpolate the data set (gt,g) at time t
Ext = interp1(tIN,Ex,t); % Interpolate the data set (gt,g) at time t
v0xt=interp1(tIN,v0x,t);
v0zt=interp1(tIN,v0z,t);

gamma1=sqrt(1-(v0xt.^2+v0zt.^2)/(cau.^2));


dvdt=[gamma1.*(-Ext+v(2).*Byt+Ext.*v(1).*v0xt/cau.^2);gamma1.*(-v(1).*Byt+Ext.*v(2).*v0xt/cau.^2)]; % 
end

