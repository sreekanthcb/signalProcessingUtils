clc;clear;close all

lh = spCommLibrary;
n         = 7;
k         = 4;
intrDepth = 8;
M         = 64;
modStyle  = 'QAM';

pktsize   = 64*k*log2(M);
n_pkts    = 1000;

snrVec    = [7:0.5:8.5];
bitsinErr = zeros(size(snrVec));

nBits     = pktsize*n_pkts;

for ii = 1:length(snrVec)
  bits_in   = round(rand(1,nBits))';

  for jj = 1:n_pkts
    frameBits = bits_in((jj-1)*pktsize+1:jj*pktsize);
    scrBits   = lh.scrambler(frameBits);
    encOut    = lh.linearBinaryBlockEncoder(scrBits,n,k);
    intrlvOut = lh.RowCol_Interleaver(encOut,intrDepth);
    symbols   = lh.bits2ModSymbols(intrlvOut,M,modStyle);

    %% Channel
    chanout   = awgn(symbols,snrVec(ii));
    %%

    decBits   = lh.modSymbols2Bits(chanout,M,modStyle);
    deintOut  = lh.RowCol_DeInterleaver(decBits,intrDepth);
    decOut    = lh.linearBinaryBlockDecoder(deintOut,n,k);
    bits_out  = lh.descrambler(decOut);

    lt = min(numel(frameBits),numel(frameBits));
    bitsinErr(ii) = bitsinErr(ii)  + sum(frameBits(:) ~= bits_out(:));
    fprintf("SNR:%d <--> Pkt id:%d (out of %d pkts) completed\n",snrVec(ii),jj,n_pkts);
  end
end

berVec = bitsinErr/nBits;

semilogy(snrVec,berVec);
