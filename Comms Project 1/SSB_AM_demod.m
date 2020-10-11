function [demodSSB] = SSB_AM_demod(m,f_cutoff,wc,t,fs, Ac)

y = m .* (Ac  * cos(wc * t));

y_low = lowpass(y,f_cutoff,fs);

demodSSB = downsample(y_low,40);

end

