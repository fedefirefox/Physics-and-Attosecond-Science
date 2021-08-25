% In this code, I am exploiting the solution of a wave equation in
% homogenous medium with second order of dispersion (see HIL book or literature or other good ultrafast-optics books)
% The units used in this code are arbitrary, however they can be easily
% adjusted to atomic units or SI units (see SPM code)


v = VideoWriter('DISPERSION.avi');
open(v);
tp=10;

vg=2 % Group velocity
kl=0.25; %k-vector (make k=w/c to use SI units)
wl=100;  %frequency arb. units
for t=-10:0.5:600 %time passed
    GVD=1; %Arb units
 z=-500:0.01:1800; %fixed space grid (greater than the plotted area)
A=tp*(sqrt(tp^2-1i*GVD*z)).^(-1).*exp(-(t-z/vg).^2./(2*(tp^2-1i*GVD*z))); %Complex envelope
E=real(A.*exp(+1i*kl*z-1i*wl*t)); %Electric Field
subplot(2,1,1)
graph1=plot(z,E);
title('Dispersion in Positive GVD Medium');
set(graph1,'LineWidth',1.5);
xlabel('Space (arb)');
ylabel('Amplitude (arb)');
axis([(-500+vg*t),(500+vg*t), -1, 1]) %Moving axis


GVD=-1;
 z=-500:0.01:1800;
A=tp*(sqrt(tp^2-1i*GVD*z)).^(-1).*exp(-(t-z/vg).^2./(2*(tp^2-1i*GVD*z)));
E=real(A.*exp(+1i*kl*z-1i*wl*t));
subplot(2,1,2)
graph2=plot(z,E);
title('Dispersion in Negative GVD Medium');
set(graph2,'LineWidth',1.5);
xlabel('Space (arb)');
ylabel('Amplitude (arb)');
axis([(-500+vg*t),(500+vg*t), -1, 1])

frame = getframe(gcf);
writeVideo(v,frame);
pause(0.05)
end
close(v)


%% NORMAL DISPERSION 
v = VideoWriter('DISPERSION2.avi');
open(v);
tp=20;
GVD=-1;

vg=2
kl=0.5;
wl=100;
for t=-10:0.5:600
    
 z=-500:0.01:1800;
A=tp*(sqrt(tp^2-1i*GVD*z)).^(-1).*exp(-(t-z/vg).^2./(2*(tp^2-1i*GVD*z)));
E=real(A.*exp(+1i*kl*z-1i*wl*t));
graph1=plot(z,E);
title('Dispersion in Positive GVD Medium');
set(graph1,'LineWidth',1.5);
xlabel('Space (arb)');
ylabel('Amplitude (arb)');
axis([(-500+vg*t),(500+vg*t), -1, 1])

frame = getframe(gcf);
writeVideo(v,frame);
end
close(v)