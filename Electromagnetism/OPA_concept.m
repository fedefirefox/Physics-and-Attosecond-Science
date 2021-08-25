% This code is particularly meaningful when put aside to the OPA discussion done in the book

% This is not a simulation of the OPA phenomenon, it is rather a conceptual
% proof of how OPA works!! 
% The key feature here is the following:
% To gain proper amplification the idler the signal and the pump must be
% spatio-temporally overlapped. The better the overlap is the greater the
% amplification. 
% The amplification stops when the PUMP is no more overlapped.


% Note that once the OPA equation have been written and solved for an easy pulse (e.g. gaussian pulses, see our book or OPA
% articles in literature) it is not that hard to simulate the process.

clear all
% Definition of group velociy for signal IDLER and PUMP, k-vector and
% frequency. Arbitrary units are used

vg1=1; 
vg2=1.5;
vg3=1;
k1=20;
k2=40;
w2=3;
w1=1.5;
w3=1.5;
k3=20;

z=-6:0.01:6; %Space where the pulse propagates 
t=0; %initial time
A=1; %initial amplitude of the 

time=-6:0.02:5;
amplitude=zeros(length(time),1);
amplitude(1)=A;
index=0;
v = VideoWriter('OPApulses.avi');
open(v);

for t=-6:0.02:5
E2=10*exp(1i*k2*z-1i*w2*t).*exp(-(z-vg2*t).^2); %PUMP PULSE
E3=1*exp(1i*k3*z-1i*w3*t).*exp(-(z-vg3*t).^2); % Signal position

subplot(3,1,1)
g1=plot(z,real(E2))
hold on
g2=plot(z,real(E3))
hold off
str = sprintf('Ovelap of Signal and Pump pulses');
title(str)
set(g1,'LineWidth',1);
set(g2,'LineWidth',1);
xlabel("Space(arb)")
    ylabel("Amplitude(arb)")
legend('PUMP','Signal')

subplot(3,1,2)
E1=A.*exp(1i*k1*z-1i*w1*t).*exp(-(z-vg1*t).^2);
g3=plot(z,real(E1),'g')
str = sprintf('Signal pulses');
title(str)
set(g3,'LineWidth',1);
xlabel("Space(arb)")
    ylabel("Amplitude(arb)")
A0=A;
A=abs(max(E2.*E3))+A0; %If the pump and the idler are overlapped A is higher. (A0 you sum the pre-existent signal)

E1=A.*exp(1i*k1*z-1i*w1*t).*exp(-(z-vg1*t).^2); %Final Signal which propagates


subplot(3,1,3)
index=index+1;
amplitude(index)=A; % You plot for this simulation time the cumulative amplification so far.
g4=plot(time+6,amplitude);
str = sprintf('Cumulative amplification in time');
title(str)
set(g3,'LineWidth',1);
xlabel("Time(arb)")
    ylabel("Amplitude(arb)")

pause(0.2);

frame = getframe(gcf);
writeVideo(v,frame);
end
close(v)