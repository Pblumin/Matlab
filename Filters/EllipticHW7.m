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

%% Elliptic

Fs = 100000;  % Sampling Frequency

Fstop1 = 15000;   % First Stopband Frequency
Fpass1 = 20000;   % First Passband Frequency
Fpass2 = 30000;   % Second Passband Frequency
Fstop2 = 35000;   % Second Stopband Frequency
Astop1 = 50;      % First Stopband Attenuation (dB)
Apass  = 5;       % Passband Ripple (dB)
Astop2 = 50;      % Second Stopband Attenuation (dB)
match  = 'both';  % Band to match exactly

h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);
Elliptic = design(h, 'ellip', 'MatchExactly', match);

n = 1024;

[H,f] = freqz(Elliptic,n,Fs);

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

sgtitle('Elliptic Magnitude-Phase Plots');

%% Applying Filter

y = filter(Elliptic,Signal);

S = fft(y);
S = fftshift(abs(S))/n;
F = Fs.*(-n/2:n/2-1)/n;

figure;
plot(F,S);
title('Fourier Transform Magnitude');
xlabel('Frequency (Hz)');
