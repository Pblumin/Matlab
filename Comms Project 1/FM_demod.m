function [demodfm] = FM_demod(m, fs, f_cutoff, w)

y = real(ifft(( (fft(m)/fs)) .* (1j * w)));

y(y < 0) = 0;  

y_low = lowpass(y,f_cutoff,fs);

demodfm = downsample(y_low,40);

end

