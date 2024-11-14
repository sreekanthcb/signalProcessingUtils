%%=========================================================================
% Synopsis     : Demonstrating Shannon Hartley Equation
% Last updated : 2024-11-07
%%=========================================================================
clc;clear;close all;

% Constant bandwidth scenerio
B  = 1;
N0 = 1;
P  = 0:10^4;
C  = B.*log2(1+P./(N0*B));
subplot(2,1,1)
plot(P,C); grid on
xlabel('power, P');ylabel('capacity, C bit/sec');
title('Capacity vs Power (Constant Bandwidth)')

% constant Power scenerio
P   = 1;
N0  = 1;
B   = 1:10^3;
C   = B.*log2(1+P./(N0*B));
subplot(2,1,2)
plot(B,C); grid on
xlabel('bandwidth, B Hz'); ylabel('capacity, C bit/sec');
title('Capacity vs Bandwidth (Constant power)')
