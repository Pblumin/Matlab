%Philip Blumin, Tony Belladonna, Paul Cucchiara
%The Cooper Union
%DSP Project 1

%% Clear Stage
clc; 
clear; 
close all;

% THE CODE BELOW IS MENT TO BE RUN SECTION BY SECTION IN ORDER AVOID
% GETTING DUPLICATE FILTER GRAPHS
%% Delta Function


%---------delta function------------

z = [1 zeros(1,3000)];
y = srconvert(z');

Y = abs(fftshift(fft(y)));
figure;
plot(Y);
title('The Low Pass Filter');
 
figure;
verify(y);
 
%% Audio file

%------------- audio file ----------------
%reading the audio file
[x,Fs] = audioread('Wagner.wav');
music = srconvert(x);
 
sound(music/2,24000);

%% Comments and write up

%the srconvert function uses the multistage method
%320 was broken up in 10,8,4

%As seen from the graphs displayed the delta function passed all the
%checkpoints for srconvert, thus the function is working properly. 

%When graphing the output of srconvert of the delta function, there is
%almost a perfect rectangle which is what we want the filter to look like

%The audio file after it is passed through srconvert sounds pretty good. 

%After the 3rd filter the number of multiplies and additions is
%signifcantly reduced, so the function is efficent

%On a system other than Matlab, such as an FPGA, the polyphase
%decomposition will be much fast