function [demodconv] = convAM_demod(conv,f_cutoff,fs)

conv(conv < 0) = 0;  

y_low = lowpass(conv,f_cutoff,fs);

demodconv = downsample(y_low,40);

end

