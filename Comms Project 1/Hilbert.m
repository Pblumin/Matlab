function output = Hilbert(input, freq)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a = -1i .* sign(freq);
input1 = fft(input);
%input1 = fftshift(input1);
o = a .* input1;

output = real(ifft(o)); 
   
end

