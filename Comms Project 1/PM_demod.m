function [demodpm] = PM_demod(Ac, m,f_cutoff,wc,t,fs)

y = Ac * m .* sin(wc * t);

y_low = lowpass(y,f_cutoff,fs);

demodpm = downsample(y_low,40);

end

