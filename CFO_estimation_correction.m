%%==================================================================================
% Synopsis     : CFO estimation and Correction Probabilistic characteristics of PAPR of a random signal
%              : Transmitted FRame is formed by along with a Preamble and midamble
% Last updated : 2021-05-11
%%==================================================================================
clc;clear;close all

nBits           = 1024;
mQam            = 16;
UpFc            = 1e6; % PAssband center frequency
bbFs            = 64e6; % Baseband Sampling rate
deltaF          = 26.81e3; % CFO introduced

%% Implementation
bitsPerConst    = log2(mQam);

bits            = round(rand(1,nBits));
bits_rearranged = reshape(bits,nBits/bitsPerConst,bitsPerConst);
M               = bi2de(bits_rearranged);
dataIQ          = qammod(M,mQam)/bitsPerConst;

% Generating Preamble
knownSeqA_Bits  = round(rand(1,256));
bits_rearranged = reshape(knownSeqA_Bits,256/2,2);
M               = bi2de(bits_rearranged);
knownSeqA_IQ    = qammod(M,4)/2; % Preamble

% Generating midamble
knownSeqB_Bits  = round(rand(1,256));
bits_rearranged = reshape(knownSeqB_Bits,256/2,2);
M               = bi2de(bits_rearranged);
knownSeqB_IQ    = qammod(M,4)/2; % Midamble

frameDataBrkP   = 212; % point at which data symbols are split to place midamble
frame           = [knownSeqA_IQ; dataIQ(1:frameDataBrkP); knownSeqB_IQ; dataIQ(frameDataBrkP+1:end)]; % Transmitted Frame

%% Channel
tDelta          = (length(knownSeqA_IQ)+frameDataBrkP)/bbFs;
max_cfoRange    = 1/(2*tDelta);
fprintf('Max possible CFO = %0.3f\nApplied CFO = %0.3f\n',max_cfoRange,deltaF);

passBandFrame   = frame .* exp(2i*pi*([UpFc+deltaF]/bbFs)*([0:length(frame)-1]')); %Passband upconversion (CFO is introduced here)

%% Reciver
rxbbFrame       = passBandFrame .* exp(2i*pi*([-UpFc]/bbFs)*([0:length(passBandFrame)-1]')); % convering back to baseband

rxrefSeqA       = rxbbFrame(1:length(knownSeqA_IQ)); % extraction of channel affected preamble
rxrefSeqB       = rxbbFrame(length(knownSeqA_IQ)+frameDataBrkP+1:length(knownSeqA_IQ)+frameDataBrkP+1+length(knownSeqA_IQ)-1);% extraction of channel affected midamble

%% CFO estimation by correlation
seqACorr        = sum(knownSeqA_IQ.*conj(rxrefSeqA));
seqBCorr        = sum(knownSeqB_IQ.*conj(rxrefSeqB));

PhaseCorrA      = angle(seqACorr);
PhaseCorrB      = angle(seqBCorr);
cfo             = (PhaseCorrA-PhaseCorrB)/(2*pi*tDelta);
fprintf('Estimated CFO = %0.3f\n',cfo);

%% CFO correction
CorrectedFrame  = rxbbFrame.* exp(2i*pi*([-cfo]/bbFs)*([0:length(passBandFrame)-1]'));

%% Plots
figure();
plot(real(frame),'-r*');hold on;plot(real(rxbbFrame),'-bo');hold on;plot(real(CorrectedFrame),'-gd');
legend('Tx Frame','Corrupted Frame','Correected Frame')
