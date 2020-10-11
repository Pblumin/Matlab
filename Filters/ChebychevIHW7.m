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
%% Chebychev I
n = 1024;
Fs = 100000;  % Sampling Frequency

Fstop = 15000;       % Stopband Frequency
Fpass = 35000;       % Passband Frequency
Astop = 40;          % Stopband Attenuation (dB)
Apass = 2;           % Passband Ripple (dB)
match = 'passband';  % Band to match exactly

h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
cheby1 = design(h, 'cheby1', 'MatchExactly', match);

[H,f] = freqz(cheby1,n,Fs);

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

sgtitle('Chebychev I Magnitude-Phase Plots');


%% Applying Filter

y = filter(cheby1,Signal);

S = fft(y);
S = fftshift(abs(S))/n;
F = Fs.*(-n/2:n/2-1)/n;

figure;
plot(F,S);
title('Fourier Transform Magnitude');
xlabel('Frequency (Hz)');

