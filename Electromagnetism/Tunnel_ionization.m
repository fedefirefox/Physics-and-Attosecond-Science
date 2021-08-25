% The following 1D and 2D code are rather naive and have been made with the
% only purpouse of showing the functional behaviour of a bendoing
% potential. 
% The real physics behind is not so simple, I suggest Prof. CD Lin book or
% Prof. Zenghu Chang On Attosecond Optics for a fresh start on Attosecond
% Science.

v = VideoWriter('Potential.avi'); % This command is used to create a video
open(v);

r=-0.10005:0.001:0.10005;
r2=-3:0.1:3;
Vatom=-(1./abs(r));
w=1;

axis manual;

for t=-2*pi:0.01:2*pi
E=1000000.*cos(w*t-r2);
Ve=-r*1000000.*cos(w*t-r);

subplot(2,1,1);

plot(r,Ve+Vatom);
axis([-0.09 0.09 -10000 10000])
grid on
xlabel("r (a.u)")
ylabel("V (a.u)")
title("Atomic Potential");


subplot(2,1,2); 
plot(r2,E);
axis([-3 3 -1000000 1000000])
xlabel("r (a.u)")
ylabel("E (a.u)")
title("Electric Field");

grid on
frame = getframe(gcf);
writeVideo(v,frame);
pause(0.01);
end
close(v);


%%%%%%%%%%%%%%%%%%%%%%%%


v = VideoWriter('Potential2d.avi');
open(v);
[x,y]=meshgrid(-0.000015:0.0000001:0.000015);
[x2,y2]=meshgrid(-3:0.1:3);

Vatom=-(1./sqrt((x.^2+y.^2)));
w=1;

axis manual;
for t=0:0.01:2*pi
Ve=-x*100000000000.*cos(w*t-x);
CO(:,:,3) = ones(25).*linspace(0,1,25); % blue
sl=surfl(x,y,Ve+Vatom)
sl.EdgeColor = 'none';
grid off;
xlabel("x (arb.u)")
ylabel("y (arb.u)")
zlabel("Amplitude (arb.u)")

title("Atomic Potential");
frame = getframe(gcf);
writeVideo(v,frame);
pause(0.01);
t
end
close(v);