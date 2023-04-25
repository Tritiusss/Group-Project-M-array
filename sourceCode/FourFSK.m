% FSK modulation
clc;
clear all;
close all;

Tb = 1;
fc1 = 2;
fc2 = 5;

t = 0:(Tb/100):Tb;

c1 = sqrt(2/Tb)*sin(2*pi*fc1*t);
c2 = sqrt(2/Tb)*sin(2*pi*fc2*t);

N = 8;
m = rand(1, N);
t1 = 0;
t2 = 1;
for i = 1:N
    t = t1:0.01:t2;
    if m(i) > 0.75
        m(i) = 01;
        ms_1 = zeros(1, length(t));
        ms_2 = ones(1, length(t));
    elseif 0.5< m(i) && m(i) <= 0.75
        m(i) = 00;
        ms_1 = zeros(1, length(t));
        ms_2 = zeros(1, length(t));
    elseif 0.25 < m(i) && m(i) <= 0.5
        m(i) = 10;
        ms_1 = ones(1, length(t));
        ms_2 = zeros(1, length(t));
    else
        m(i) = 11;
        ms_1 = ones(1, length(t));
        ms_2 = ones(1, length(t));
    end
    fsk_sig1(i,:) = c1.*ms_1;
    fsk_sig2(i,:) = c2.*ms_2;
    fsk=fsk_sig1+fsk_sig2;
    subplot(3, 2, 4);
    plot(t, fsk(i,:));
    xlabel('time(sec)');
    ylabel('s(t)');
    grid on;
    hold on;
    t1 = t1 + (Tb+.01);
    t2 = t2 + (Tb+.01);
end
hold off
subplot(3, 2, 1);
stem(m);
title('message data');
xlabel('n---->');
ylabel('b(n)');
grid on;
subplot(3, 2, 2);
plot(t, c1);
title('carrier signal-1');
xlabel('time(sec)');
ylabel('cl(t)');
grid on;
subplot(3, 2, 3);
plot(t, c2);
xlabel('time(sec)');
ylabel('cl(t)');
grid on;


%both message symbol 11 and 00 will make x equal to 0. Therefore I cannot divided
%them into two parts. And 0 and 1 are correct area
%demodulation
t3 = 0;
t4 = Tb;
for i = 1:N
    t=t3:(Tb/100):t4;
    x1 = sum(c1.*fsk_sig1(i,:));
    x2 = sum(c2.*fsk_sig2(i,:));
    x = x1 - x2;
    if x >= 0.5
        demod(i)= 10;
    elseif x == 0
        demod(i)= 11;
    elseif x <= -0.5
        demod(i)= 1;
    else
        demod(i)=0;
    end
    t3 = t3+(Tb+.01);
    t4 = t4+(Tb+.01);
end
subplot(3, 2, 5);
stem(demod);
title('demodulated data');
xlabel('n---->');
ylabel('b(n)');
grid on;


disp(m)