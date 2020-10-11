% Philip Blumin
% MATLab Section B HW 7

%% Clear Stage
clc;
clear;
close all;

%% Filter Design
filterDesigner

%% Signal
t = linspace(0,2,1024);
f = reshape(1:50000,50000,1);
Signal = sum(sin(2 * pi * f .* t));
%% BytterWorth
n = 1024;
Fs = 100000;  % Sampling Frequency

Fpass = 10000;       % Passband Frequency
Fstop = 20000;       % Stopband Frequency
Apass = 5;           % Passband Ripple (dB)
Astop = 50;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly


h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
ButterWorth = design(h, 'butter', 'MatchExactly', match);

[H,f] = freqz(ButterWorth,n,Fs);

hdb = 20*(log10(abs(H)));
hph = unwrap(angle(H))* 180/pi;

figure;
subplot(2,1,1);
plot(f,hdb);
title('Magnitude Response in DB');
xlabel('w (Hz)');

subplot(2,1,2);
plot(f,hph);
title('Phase Response');
xlabel('w (Hz)');

sgtitle('Butter Worth Magnitude-Phase Plots');

%% Applying Filter

y = filter(ButterWorth,Signal);

S = fft(y);
S = fftshift(abs(S))/n;
F = Fs.*(-n/2:n/2-1)/n;

figure;
plot(F,S);
title('Fourier Transform Magnitude');
xlabel('Frequency (Hz)');


