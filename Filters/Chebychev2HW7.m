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

%% Chebychev II
Fs = 100000;  % Sampling Frequency

Fpass1 = 5000;        % First Passband Frequency
Fstop1 = 15000;       % First Stopband Frequency
Fstop2 = 35000;       % Second Stopband Frequency
Fpass2 = 45000;       % Second Passband Frequency
Apass1 = 5;           % First Passband Ripple (dB)
Astop  = 50;          % Stopband Attenuation (dB)
Apass2 = 5;           % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly

h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, Apass2, Fs);
cheby2 = design(h, 'cheby2', 'MatchExactly', match);

n = 1024;
[H,f] = freqz(cheby2,n,Fs);

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

sgtitle('Chebychev II Magnitude-Phase Plots');

%% Applying Filter

y = filter(cheby2,Signal);

S = fft(y);
S = fftshift(abs(S))/n;
F = Fs.*(-n/2:n/2-1)/n;

figure;
plot(F,S);
title('Fourier Transform Magnitude');
xlabel('Frequency (Hz)');
