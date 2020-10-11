function [u_ussb] = SSB_AM_modulate(mt,a,Ac,f,fs,wc)

N = length(mt);
T = N/fs;
N2 = N * 40; %upsampling N
t = linspace(0,T,N);
t2 = linspace(0,T,N2); %upsampled t

mt_scaled = mt ./ a;

mt_upscaled = interp1(t, mt_scaled, t2);

m_hat = Hilbert(mt_upscaled, f);

u_ussb = (Ac * mt_upscaled .* cos(wc * t2)) + (Ac * m_hat .* sin(wc * t2));

u_ussb = u_ussb + cos(wc * t2);

end