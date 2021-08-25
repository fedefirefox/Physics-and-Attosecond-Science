% This is another simple code, where I am just exploting the concept of Fourier Transform to show why it is impossible to propagate a "pulse" % which duration is below the optical cycle. 


v = VideoWriter('Single_cycle.avi');
open(v);

for sigma=1:-0.0005:0.2
freq=5;

Fs = 100;           % Sampling frequency
t= -3:1/Fs:3;  % Time vector 
L = length(t);      % Signal length
n = 2^nextpow2(L);

E3=exp(1i*2*pi*freq*t+1i*pi).*exp(-(t/sigma).^2/2);

subplot(2,1,1)
g1=plot(t,real(E3))
    set(g1,'LineWidth',1);
    xlabel("Time(fs)")
    ylabel("Amplitude(arb)")
    str = sprintf('Electric FIELD \\lambda= 60nm time domain t_p= %f fs', sigma);
    title(str);
    

Y = fft(real(E3),n);
f = Fs*(0:(n/2))/n;
P = abs(Y/n);

subplot(2,1,2)
g2=plot(f,P(1:n/2+1),'r')

    set(g2,'LineWidth',1);
    xlabel("Frequency(PHz)")
    ylabel("Amplitude(arb)")
    str = sprintf('Electric FIELD, frequency domain');
    title(str);


frame = getframe(gcf);
writeVideo(v,frame);
end

for sigma=0.2:-0.0005:0.01
freq=5;

Fs = 100;           % Sampling frequency
t= -3:1/Fs:3;  % Time vector 
L = length(t);      % Signal length
n = 2^nextpow2(L);

E3=exp(1i*2*pi*freq*t+1i*pi).*exp(-(t/sigma).^2/2);

subplot(2,1,1)
g1=plot(t,real(E3))
    set(g1,'LineWidth',1);
    xlabel("Time(fs)")
    ylabel("Amplitude(arb)")
    str = sprintf('Electric FIELD \\lambda= 60nm time domain t_p= %f fs', sigma);
    title(str);
    

Y = fft(real(E3),n);
f = Fs*(0:(n/2))/n;
P = abs(Y/n);

subplot(2,1,2)
g2=plot(f,P(1:n/2+1),'r')

    set(g2,'LineWidth',1);
    xlabel("Frequency(PHz)")
    ylabel("Amplitude(arb)")
        str = sprintf('Electric FIELD \\lambda= 60nm beyond optical cycle');
    title(str);
frame = getframe(gcf);
writeVideo(v,frame);

end
close(v);
