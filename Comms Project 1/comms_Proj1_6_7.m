% Philip Blumin
% Comms Project 1
% with noise

%% Clear Stage
clc;
clear;
close all;

%% Scaling Message
%WARNING WHEN LISTENING TO MP3 FILE PLZ LOWER VOLUME FOR SAFTY PURPOSES
%ALL MAGNITUDE PLOTS HAVE Y AXIS IN LOG SCALE USING SEMILOGYY

[mt, fs] = audioread('Carl_Smith_Comms.mp3'); 
mt = mt(:, 1)';

N = length(mt);
T = N/fs;
N2 = N * 40; %upsampling N
fs2 = fs * 40; %upsampling fs
t = linspace(0,T,N);
t2 = linspace(0,T,N2); %upsampled t
max = max(abs(mt));
wd = linspace(-pi,pi,N); %wd
wd2 = linspace(-pi,pi,N2); %wd
f = (wd * fs)/(2 * pi);
f2 = (wd2 * fs2)/(2 * pi);
Ac = 1;

mt_scaled = mt ./ max;

mt_upscaled = interp1(t, mt_scaled, t2);

MT_upscaled = fft(mt_upscaled);

MT_upscaled = fftshift(MT_upscaled/fs2);

fc = 600E3; %carrier freq. I decided on 600khz after googling carrier freqs used for AM

wc = 2 * pi * fc;

f_cutoff = fs/2;
variance = 0.05 ;

%modulation index multiplier
a = [2 1 0.5];

figure;
plot(t2,mt_upscaled);
xlabel('t (sec)');
title('Original Signal');
ylabel('Magnitude');

%as seen the message is now scaled so that the maximum abs is 1


%% Convenional AM

%modulate
for i = 1:3
    index = a(i);
    [conv] = convAM_modulate(mt,max * a(i),Ac,fs,wc);

    CONV = fft(conv);

    CONV = fftshift(CONV/fs2);

    %adding noise

    %[conv_noise] = convAM_modulate(mt_noise,max,Ac,fs,wc);
    noisevec = sqrt(variance)*randn(length(mt_upscaled),1); %white noise
    noisevec = noisevec'; 

    r_conv = snr(conv,noisevec);

    disp("SNR for conv with modulation index factor of " + index + ": " + r_conv);

    conv_noise = conv + noisevec;

    CONV_noise = fft(conv_noise);

    CONV_noise = fftshift(CONV_noise/fs2);

    %demod
    [demodconv] = convAM_demod(conv_noise,f_cutoff,fs2);

    demodconv = lowpass(demodconv,5000,fs2);

    Demodconv = fft(demodconv);

    Demodconv = fftshift(Demodconv/fs2);
    
    figure;
    subplot(2,2,1);
    semilogy(f2,abs(MT_upscaled));
    xlabel('w (Hz)');
    title('Original Signal');
    ylabel('Magnitude (log)');

    subplot(2,2,2);
    semilogy(f2, abs(CONV));
    title('AM modulated');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,3);
    semilogy(f2, abs(CONV_noise));
    title('AM modulated with noise');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,4);
    semilogy(f, abs(Demodconv));
    title('AM demodulated');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');
    
    titleconv = ['Modulation index factor = ',num2str(index)];
    sgtitle(titleconv);
    
end
% figure;
% plot(t,demodconv);
% title('conventional AM demodulated in time domain');
% xlabel('t (sec)');

%sound(demodconv, fs);

%% PM

for i = 1:3
    index = a(i);
    %modulate
    [pm] = PM_mod(mt,Ac, wc, fs, index * 1);

    PM = fft(pm);

    PM = fftshift(PM/fs2);

    %2 deltas with just regular plot in freq domain

    %adding noise

    noisevec = sqrt(variance)*randn(length(mt_upscaled),1); %white noise
    noisevec = noisevec'; 
    
    r_pm = snr(pm,noisevec);
    
    disp("SNR for pm with modulation index factor of " + index + ": " + r_pm);
    
    pm_noise = pm + noisevec;

    PM_noise = fft(pm_noise);

    PM_noise = fftshift(PM_noise/fs2);

    %demodulation

    [demodpm] = PM_demod(Ac, pm_noise ,f_cutoff,wc,t2,fs2);

    demodpm = lowpass(demodpm,5000,fs2);

    demodPM = fft(demodpm);

    demodPM = fftshift(demodPM/fs2);

    figure;
    subplot(2,2,1);
    semilogy(f2,abs(MT_upscaled));
    xlabel('w (Hz)');
    title('Original Signal');
    ylabel('Magnitude (log)');
    
    subplot(2,2,2);
    semilogy(f2, abs(PM));
    title('PM modulated in frequency domain');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,3);
    semilogy(f2, abs(PM_noise));
    title('PM modulated with noise');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,4);
    semilogy(f, abs(demodPM));
    title('PM demodulated in frequency domain');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');
    
    titlepm = ['Modulation index factor = ',num2str(index)];
    sgtitle(titlepm);

end

% figure;
% plot(t,demodpm);
% title('PM demodulated in time domain');
% xlabel('t (sec)');

%soundsc(demodpm, fs);

%% FM
Ac1 = 15;
%modulate
wc2 = 2 * pi * 1000E3;

for i = 1:3
    index = a(i);

    [fm] = FM_mod(mt,Ac, wc2, fs, index * 100000);

    FM = fft(fm);

    FM = fftshift(FM/fs2);

    %adding noise

    noisevec = sqrt(variance)*randn(length(mt_upscaled),1); %white noise
    noisevec = noisevec'; 
    
    r_fm = snr(fm,noisevec);
    
    disp("SNR for fm with modulation index factor of " + index + ": " + r_fm);
    
    fm_noise = fm + noisevec;

    FM_noise = fft(fm_noise);

    FM_noise = fftshift(FM_noise/fs2);

    %demod

    w = f2 * (2 * pi);

    [demodfm] = FM_demod(fm_noise, fs2, f_cutoff, w);

    demodfm = lowpass(demodfm,5000,fs2);

    demodFM = fft(demodfm);

    demodFM = fftshift(demodFM/fs2);

    figure;
    subplot(2,2,1);
    semilogy(f2,abs(MT_upscaled));
    xlabel('w (Hz)');
    title('Original Signal');
    ylabel('Magnitude (log)');
    
    subplot(2,2,2);
    semilogy(f2, abs(FM));
    title('FM modulated in frequency domain');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,3);
    semilogy(f2, abs(FM_noise));
    title('FM modulated with noise');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,4);
    semilogy(f, abs(demodFM));
    title('FM demodulated in frequency domain');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');
    
    titlefm = ['Modulation index factor = ',num2str(index)];
    sgtitle(titlefm);

end
% 
% figure;
% plot(t,demodfm);
% title('FM demodulated in time domain');
% xlabel('t (sec)');

%soundsc(demodfm, fs);
