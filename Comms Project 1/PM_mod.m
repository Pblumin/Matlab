function [pm] = PM_mod(mt,Ac, wc, fs, k)

N = length(mt);
T = N/fs;
N2 = N * 40; %upsampling N
t = linspace(0,T,N);
t2 = linspace(0,T,N2); %upsampled t

a = max(abs(mt));

mt_scaled = mt ./ a;

mt_upscaled = interp1(t, mt_scaled, t2);

phase = k * mt_upscaled;

pm = Ac * cos(wc * t2 + phase);

end

