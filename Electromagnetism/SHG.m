% In the simulation, all points of the space participate in the "generation" of a sinusoidal second harmonic wave,
% with a weight that depends on the local value of the non linear polarization of the medium at 2ω. 

% The simulation shows that a nice
% propagative mode is obtained in the phase-matching condition.
% ONLY in this condition, all the points in space generate a wave that "sums
% coherently".

% However, it does not clearly show
% why the other modes do not propagate, due to the approximations used.

% As always a proper resolution of Maxwell equation would provide a clear
% answer on the phenomenon


clear all
v = VideoWriter('PHASE_MATCHING.avi');
open(v);
z=0:0.1:100;
t=0;
w1=0.8;
k1=1;
k2=5;
P1=sin(k1*z-w1*t);
zgen=0:0.1:100;
k2v=[1,1.5,2,2,2.5,3,4,5,6,7];
for j=1:length(k2v)
for t=0:0.1:10
    k2=k2v(j)
    Etot=0;
P1=sin(2*k1*z-2*w1*t);
E=sin(k1*z-w1*t);
 for i=1:length(zgen)
     Ea(i)=P1(1,i)*sin(k2*(z(i))-2*w1*t);
 end
 subplot(2,1,1) 
 g1=plot(z,E);
 str = sprintf('Electric FIELD, INITIAL FIELD k_{w} = 1', k2);
title(str)
set(g1,'LineWidth',1);
xlabel("Space(arb)")
    ylabel("Amplitude(arb)")
    
  subplot(2,1,2)
    g2=plot(z,Ea);
    set(g2,'LineWidth',1);
    xlabel("Space(arb)")
    ylabel("Amplitude(arb)")
    str = sprintf('Electric FIELD, Second HARMONIC k_{2w} = %f', k2);
    if k2==2
        str=sprintf('k_{2w} is %f PHASE MATCHING', k2);
    end
title(str)
       pause(0.1);
       frame = getframe(gcf);
writeVideo(v,frame);
end

end
close(v)