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

variance = [0.25 0.1 0.05];

figure;
plot(t2,mt_upscaled);
xlabel('t (sec)');
ylabel('Magnitude');
title('Original Signal');

%as seen the message is now scaled so that the maximum abs is 1


%% SSB AM

%modulate
[u_ussb] = SSB_AM_modulate(mt,max,Ac,f2,fs,wc);

U_USSB = fft(u_ussb);

U_USSB = fftshift(U_USSB/fs2);

%adding noise
%[u_ussb_noise] = SSB_AM_modulate(mt_noise,max,Ac,f2,fs,wc);
for i = 1:3
    a = variance(i);
    noisevec = sqrt(variance(i))*randn(length(mt_upscaled),1); %white noise
    noisevec = noisevec'; 
    
    r_ssb = snr(u_ussb,noisevec);
    
    disp("SNR for USSB (var = " + a + "): " + r_ssb);
    
    u_ussb_noise = u_ussb + noisevec;
    
    U_ussb_noise = fft(u_ussb_noise);

    U_ussb_noise = fftshift(U_ussb_noise/fs2);
    

    %demod
    [demodSSB] = SSB_AM_demod(u_ussb_noise,f_cutoff,wc,t2,fs2,Ac);

    demodSSB = lowpass(demodSSB,5000,fs2);

    DemodSSB = fft(demodSSB);

    DemodSSB = fftshift(DemodSSB/fs2);

    figure;
    subplot(2,2,1);
    semilogy(f2,abs(MT_upscaled));
    xlabel('w (Hz)');
    title('Original Signal')
    ylabel('Magnitude (log)');
    
    subplot(2,2,2);
    semilogy(f2, abs(U_USSB));
    title('SSB modulated');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,3);
    semilogy(f2, abs(U_ussb_noise));
    title('SSB modulated with noise');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,4);
    semilogy(f, abs(DemodSSB));
    title('SSB demodulated');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');
    
    titlessb = ['variance = ',num2str(a)];
    sgtitle(titlessb);
end

% figure;
% plot(t,demodSSB);
% title('upper SSB demodulated in time domain');
% xlabel('t (sec)');

soundsc(demodSSB, fs);

%% Convenional AM

%modulate
[conv] = convAM_modulate(mt,max,Ac,fs,wc);

CONV = fft(conv);

CONV = fftshift(CONV/fs2);

%adding noise

%[conv_noise] = convAM_modulate(mt_noise,max,Ac,fs,wc);
for i = 1:3
    a = variance(i);
    noisevec = sqrt(variance(i))*randn(length(mt_upscaled),1); %white noise
    noisevec = noisevec'; 
    
    r_conv = snr(conv,noisevec);
    
    disp("SNR for conv (var = " + a + "): " + r_conv);
    
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
    title('Original Signal')
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
    
    titleconv = ['variance = ',num2str(a)];
    sgtitle(titleconv);
end

% figure;
% plot(t,demodconv);
% title('conventional AM demodulated in time domain');
% xlabel('t (sec)');

sound(demodconv, fs);

%% PM

%modulate
[pm] = PM_mod(mt,Ac, wc, fs, 1);

PM = fft(pm);

PM = fftshift(PM/fs2);

%2 deltas with just regular plot in freq domain

%adding noise

for i = 1:3
    a = variance(i);
    noisevec = sqrt(variance(i))*randn(length(mt_upscaled),1); %white noise
    noisevec = noisevec'; 
    
    r_pm = snr(pm,noisevec);
    
    disp("SNR for pm (var = " + a + "): " + r_pm);
    
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
    title('PM modulated');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,3);
    semilogy(f2, abs(PM_noise));
    title('PM modulated with noise');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,4);
    semilogy(f, abs(demodPM));
    title('PM demodulated');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');
    
    titlepm = ['variance = ',num2str(a)];
    sgtitle(titlepm);

end

% figure;
% plot(t,demodpm);
% title('PM demodulated in time domain');
% xlabel('t (sec)');

soundsc(demodpm, fs);

%% FM
Ac1 = 15;
%modulate
wc2 = 2 * pi * 1000E3;
[fm] = FM_mod(mt,Ac, wc2, fs, 100000);

FM = fft(fm);

FM = fftshift(FM/fs2);

%adding noise

%[fm_noise] = FM_mod(mt_noise,Ac, wc2, fs, max, 100000);
for i = 1:3
    a = variance(i);
    noisevec = sqrt(variance(i))*randn(length(mt_upscaled),1); %white noise
    noisevec = noisevec'; 
    
    r_fm = snr(fm,noisevec);
    
    disp("SNR for fm (var = " + a + "): " + r_fm);
    
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
    title('FM modulated');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,3);
    semilogy(f2, abs(FM_noise));
    title('FM modulated with noise');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');

    subplot(2,2,4);
    semilogy(f, abs(demodFM));
    title('FM demodulated');
    xlabel('w (Hz)');
    ylabel('Magnitude (log)');
    
    titlepm = ['variance = ',num2str(a)];
    sgtitle(titlepm);

end
% 
% figure;
% plot(t,demodfm);
% title('FM demodulated in time domain');
% xlabel('t (sec)');

soundsc(demodfm, fs);
