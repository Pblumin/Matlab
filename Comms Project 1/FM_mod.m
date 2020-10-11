function [fm] = FM_mod(mt,Ac, wc, fs, k)
%kf

N = length(mt);
T = N/fs;
N2 = N * 40; %upsampling N
fs2 = fs * 40; %upsampling fs
t = linspace(0,T,N);
t2 = linspace(0,T,N2); %upsampled t

a = max(abs(mt));

mt_scaled = mt ./ a;

mt_upscaled = interp1(t, mt_scaled, t2);

phase = 2 * pi * k * (cumsum(mt_upscaled)/(fs2));

fm = Ac * cos(wc * t2 + phase); 

end

