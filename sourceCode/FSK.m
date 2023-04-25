% FSK modulation
clc;
clear all;
close all;

Tb = 1;
fc1 = 2;
fc2 = 5;

t = 0:(Tb/100):Tb;

% Generate carrier signal
c1 = sqrt(2/Tb)*sin(2*pi*fc1*t);
%c1n = c1 + sqrt(t).*randn(size(t)) + 1;
c1n = awgn(c1, 10, 'measured');
c2 = sqrt(2/Tb)*sin(2*pi*fc2*t);
c2n = awgn(c2, 10, 'measured');

% Generate random message signal
N = 8;
m = rand(1, N);
t1 = 0;
t2 = Tb;
for i = 1:N
    t = t1:(Tb/100):t2;
    if m(i)>0.5
        m(i) = 1;
        m_s1 = ones(1, length(t));
        m_s2 = zeros(1, length(t));
    else
        m(i) = 0;
        m_s1 = zeros(1, length(t));
        m_s2 = ones(1, length(t));
    end
    message(i,:) = m_s1;
    fsk_sig1(i,:) = c1n.*m_s1;
    fsk_sig2(i,:) = c2n.*m_s2;
    fsk = fsk_sig1+fsk_sig2;
    subplot(3, 2, 4);
    plot(t, fsk(i,:));
    title('FSK signal');
    xlabel('time(sec)');
    ylabel('s(t)');
    grid on;
    hold on;
    t1 = t1 + (Tb+.01);
    t2 = t2 + (Tb+.01);
end
hold off
% plot carrier signal
subplot(3, 2, 1); 
stem(m);
title('binary data');
xlabel('n---->');
ylabel('b(n)');
grid on;
subplot(3, 2, 2);
plot(t, c1n);
title('noise carrier signal-1');
xlabel('time(sec)');
ylabel('cl(t)');
grid on;
subplot(3, 2, 3);
plot(t, c2n);
title('noise carrier signal-2');
xlabel('time(sec)');
ylabel('c2(t)');
grid on;

%demodulation
t3 = 0;
t4 = Tb;
for i=1:N
    t=t3:Tb/100:t4;
    x1 = sum(c1.*fsk_sig1(i,:));
    x2 = sum(c2.*fsk_sig2(i,:));
    x = x1 - x2;
    if x > 0
        demod(i) = 1;
    else
        demod(i) = 0;
    end
    t3 = t3+(Tb+.01);
    t4 = t4+(Tb+.01);
end
subplot(3, 2, 5);
stem(demod);
title('demodulated data');
xlabel('n');
ylabel('b(n)');
grid on;
disp(m)
