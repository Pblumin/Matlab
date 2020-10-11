function [conv] = convAM_modulate(mt,a,Ac,fs,wc)
%a is the modulation index; k = 1 so mod index will be just max

N = length(mt);
T = N/fs;
N2 = N * 40; %upsampling N
t = linspace(0,T,N);
t2 = linspace(0,T,N2); %upsampled t

mt_scaled = mt ./ a;

mt_upscaled = interp1(t, mt_scaled, t2);

c = cos(wc * t2);

conv = Ac * (1 + mt_upscaled) .* c;

end

