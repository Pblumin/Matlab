function [out] = srconvert(in)

%Need to upsample by 147 and downsample by 320
%Multistage split up into 10, 8, 4

M = 147;

%---------cheby filter 1----------
x_up = upsample(in,10);

Wp = 1/10;
Ws = 1.2 * Wp;
Rp = 0.1;
Rs = 75.4;

[n,Wn] = cheb2ord(Wp,Ws,Rp,Rs);
[b,a] = cheby2(n,Rs,Wn);

out1 = filter(b,a,x_up);

figure;
freqz(b,a);
title('cheby filter 1');

%Using formulas from Oppenheim, section 4.7.5, we see that the size of the
%filter is N = 18. Since we are upsampling by 10, L = 10. The amount of
%additions per unit time are calculated by N*L-1 and the multiplications
%are calculated by N*L.
disp('------------------------------------------------------------------');
disp('For the first stage, the size of the Cheby filter is N = 18');
disp('There are 179 additions and 180 multiplications per unit time')
disp('------------------------------------------------------------------');

%---------cheby filter 2----------
x_up1 = upsample(out1,8);

Wp1 = 1/8;
Ws1 = 1.2 * Wp1;
Rp1 = 0.03;
Rs1 = 73;

[n1,Wn1] = cheb2ord(Wp1,Ws1,Rp1,Rs1);
[b1,a1] = cheby2(n1,Rs1,Wn1);

out2 = filter(b1,a1,x_up1);

figure;
freqz(b1,a1);
title('cheby filter 2');

%Using formulas from Oppenheim, section 4.7.5, we see that the size of the
%filter is N = 19. Since we are upsampling by 8, L = 8. The amount of
%additions per unit time are calculated by N*L-1 and the multiplications
%are calculated by N*L.
disp('------------------------------------------------------------------');
disp('For the second stage, the size of the Cheby filter is N = 19');
disp('There are 151 additions and 152 multiplications per unit time')
disp('------------------------------------------------------------------');

%---------cheby filter 3----------
x_up2 = upsample(out2,4);

Wp2 = 1/4;
Ws2 = 1.2 * Wp2;
Rp2 = 0.1;
Rs2 = 74.193;

[n2,Wn2] = cheb2ord(Wp2,Ws2,Rp2,Rs2);
[b2,a2] = cheby2(n2,Rs2,Wn2);

[h2,t] = impz(b2,a2);

%---------polyphase----------
E = poly1(h2',4);
E = upsample(E',4)'; %ready to be offset

polyout = conv(E(1,:),x_up2);
for i = 2:4
    polyout = polyout';
    polyout = [polyout 0];
    polyout = polyout';
    polyout = polyout + conv([zeros(1,i-1), E(i,:)],x_up2);
end

figure;
freqz(b2,a2);
title('cheby3 filter 3');

%Using formulas from Oppenheim, section 4.7.5, we see that the size of the
%filter is N = 17. Since we are upsampling by 4, L = 4. The amount of
%additions per unit time are calculated by L*(N/L-1) and the multiplications
%are calculated by L*(N/L).
disp('----------------------------------------------------------------------------------------');
disp('For the third stage (Polyphase implementation), the size of the Cheby filter is N = 17');
disp('There are 13 additions and 17 multiplications per unit time')
disp('----------------------------------------------------------------------------------------');


%---------final output downsample----------
out = downsample(polyout,147);
out = out./max(abs(out));

end