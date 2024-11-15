%%==================================================================================================
% Synopsis     : Timing recovery of Frame using Auto correlation
%              : Frame is formed having preamble which is a two time repeatation of a known pattern
% Last updated : 2019-02-18
%%==================================================================================================
clc;clear;close all;

%% Input parameters
frameLen       = 1024;
PreambleLen    = 64;
timeError      = 23; % no of random samples added before frame to introduce timing error

%% Implementation
DataLen        = frameLen - PreambleLen;
n_preambleBits = PreambleLen/2; % preamble has similar two halves
n_DataBits     = DataLen;

dataBits       = round(rand(1,n_DataBits));
preambleBits   = round(rand(1,n_preambleBits));

data           = 2*dataBits-1; % BPSK modulation
preamble       = 2*preambleBits-1;

FRAME          = [preamble preamble data]; % formation of FRAME with preamble and DAta

randSamples    = 2*round(rand(1,timeError))-1;
FRAME_inError  = [randSamples FRAME]; % time corrupted FRAME

% Timing Recovery
corrOut        = xcorr(FRAME_inError);
corrOut_mag    = corrOut.*conj(corrOut); % magnitude of correlation output
[~,pos]        = max(corrOut_mag); % finding peak of corrleation output
preambleStartPoint = pos - frameLen + 1; % identifying start of FRAME

Extracted_FRAME    = FRAME_inError(preambleStartPoint:end);
Extracted_Data     = Extracted_FRAME(PreambleLen+1:end);

if(Extracted_Data == data)
  fprintf("Timing Recovery succesfull\nExtracted Data is same as the Transmitted Data Samples\n");
else
  fprintf("Timing recovery Failed\n");
end

