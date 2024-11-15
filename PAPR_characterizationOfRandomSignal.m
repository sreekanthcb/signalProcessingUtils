%%==================================================================================
% Synopsis     : Analysing Probabilistic characteristics of PAPR of a random signal
% Last updated : 2019-02-18
%%==================================================================================
clc;clear;close all

f           = 45e3; % frequency of Sinetone
fs          = 1e6; % sampling rate
n_samples   = 1e5;  % no of samples
papr_window = 0:0.1:14;
snrValue    = 1;

%% Implementation
sinetoneClean      = sin(2*pi*f/fs*(0:n_samples-1));
sinewithNoise      = awgn(sinetoneClean,snrValue);

pow_cleanSine      = sinetoneClean.*conj(sinetoneClean)/mean(sinetoneClean.*conj(sinetoneClean)); % ratio of instanteneous power of a sample to average power of signal
pow_dB_cleanSine   = 10*log10(pow_cleanSine); % above in decibel
pow_noisySine      = sinewithNoise.*conj(sinewithNoise)/mean(sinewithNoise.*conj(sinewithNoise)); % ratio of instanteneous power of a sample to average power of signal
pow_dB_noisySine   = 10*log10(pow_noisySine); % above in decibel

samplesAbovePAPR_clean = [];
samplesAbovePAPR_noisy = [];

for ii = 1:length(papr_window)
    samplesAbovePAPR_clean(ii) = sum(pow_dB_cleanSine>papr_window(ii))/n_samples; % percentage of no of samples having power above papr_window(ii)
    samplesAbovePAPR_noisy(ii) = sum(pow_dB_noisySine>papr_window(ii))/n_samples; % same for noisy sine
end

%ploting result
semilogy(papr_window,samplesAbovePAPR_clean,'-bo');grid on;hold on
semilogy(papr_window,samplesAbovePAPR_noisy,'-r*');
title('Probability of given peak power');
ylabel('Probability');
xlabel('power level above mean power');
legend('Clean Sine tone','Noisy Sine tone');
