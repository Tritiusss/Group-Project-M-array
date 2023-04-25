clc;
clear all;
close all;

Tb = 1;
fc1 = 2;
fc2 = 5;
fc3 = 8;

t = 0:(Tb/100):Tb;

c1 = sqrt(2/Tb)*sin(2*pi*fc1*t);
c2 = sqrt(2/Tb)*sin(2*pi*fc2*t);
c3 = sqrt(2/Tb)*sin(2*pi*fc3*t);

% random message signal generate
N = 8;
m = rand(1, N);
t1 = 0;
t2 = 1;
for i = 1:N
    t = t1:Tb/100:t2;
    if m(i) > 0.875
        m(i) = 000;
        ms_1 = zeros(1, length(t));
        ms_2 = zeros(1, length(t));
        ms_3 = zeros(1, length(t));
    elseif 0.75 < m(i) && m(i) <= 0.875
        m(i) = 001;
        ms_1 = zeros(1, length(t));
        ms_2 = zeros(1, length(t));
        ms_3 = ones(1, length(t));
    elseif 0.625 < m(i) && m(i) <= 0.75
        m(i) = 010;
        ms_1 = zeros(1, length(t));
        ms_2 = ones(1, length(t));
        ms_3 = zeros(1, length(t));
    elseif 0.5 < m(i) && m(i) <= 0.625
        m(i) = 100;
        ms_1 = ones(1, length(t));
        ms_2 = zeros(1, length(t));
        ms_3 = zeros(1, length(t));
    elseif 0.375 < m(i) && m(i) <= 0.5
        m(i) = 011;
        ms_1 = zeros(1, length(t));
        ms_2 = ones(1, length(t));
        ms_3 = ones(1, length(t));
    elseif 0.25 < m(i) && m(i) <= 0.375
        m(i) = 101;
        ms_1 = ones(1, length(t));
        ms_2 = zeros(1, length(t));
        ms_3 = ones(1, length(t));
    elseif 0.125 < m(i) && m(i) <= 0.25
        m(i) = 110;
        ms_1 = ones(1, length(t));
        ms_2 = ones(1, length(t));
        ms_3 = zeros(1, length(t));
    else
        m(i) = 111;
        ms_1 = ones(1, length(t));
        ms_2 = ones(1, length(t));
        ms_3 = ones(1, length(t));
    end
    fsk_sig1(i,:) = c1.*ms_1;
    fsk_sig2(i,:) = c2.*ms_2;
    fsk_sig3(i,:) = c3.*ms_3;
    fsk = fsk_sig1 + fsk_sig2 + fsk_sig3;
    subplot(3, 2, 4);
    plot(t, fsk(i, :));
    title('8-FSK signal')
    xlabel('time(sec)');
    ylabel('s(t)');
    grid on;
    hold on;
    t1 = t1 + (Tb+.01);
    t2 = t2 + (Tb+.01);
end
hold off
subplot(3, 2, 1);
plot(t, c1);
title('carrier signal-1');
xlabel('time(sec)');
ylabel('cl(t)');
grid on;
subplot(3, 2, 2);
plot(t, c2);
title('carrier signal-2');
xlabel('time(sec)');
ylabel('cl(t)');
grid on;
subplot(3, 2, 3);
plot(t, c3);
title('carrier signal-3');
xlabel('time(sec)');
ylabel('cl(t)');
grid on;
disp(m)

