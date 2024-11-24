%%==================================================================================
% Synopsis     : Signal Processing and Communication utilities Library
%              : Example usecase can be found in the accompanying testbench file
% Last updated : 2024-11-20
%%==================================================================================

function libHandle = spCommLibrary

  libHandle.bits2Msg                  = @bits2Msg;
  libHandle.msg2Bits                  = @msg2Bits;
  libHandle.RowCol_Interleaver        = @RowCol_Interleaver;
  libHandle.RowCol_DeInterleaver      = @RowCol_DeInterleaver;
  libHandle.turboLikeInterleaver      = @turboLikeInterleaver;
  libHandle.turboLikeDeInterleaver    = @turboLikeDeInterleaver;
  libHandle.scrambler                 = @scrambler;
  libHandle.descrambler               = @descrambler;
  libHandle.spectrumVisualizer        = @spectrumVisualizer;
  libHandle.PolyPhaseInterpolator     = @PolyPhaseInterpolator;
  libHandle.PolyPhaseDecimator        = @PolyPhaseDecimator;
  libHandle.findEVM                   = @findEVM;
  libHandle.bits2ModSymbols           = @bits2ModSymbols;
  libHandle.modSymbols2Bits           = @modSymbols2Bits;
  libHandle.linearBinaryBlockEncoder  = @linearBinaryBlockEncoder;
  libHandle.linearBinaryBlockDecoder  = @linearBinaryBlockDecoder;

end


function msg = bits2Msg(bits,nBitsPerMsg)

if nargin == 1
    nBitsPerMsg = 4;
end

msg = bi2de(reshape(bits,nBitsPerMsg,[])','left-msb');
msg = msg(:)';
end

function bits3 = msg2Bits(msg,nBitsPerMsg)

if nargin == 1
    nBitsPerMsg = 4;
end

bits1 = de2bi(msg,nBitsPerMsg,'left-msb');
bits2 = bits1';
bits3 = bits2(:)';

end

function out = RowCol_Interleaver(inp,k)

out = reshape(inp,[],length(inp)/k)';
out = out(:)';

end

function out = RowCol_DeInterleaver(inp,k)

out = reshape(inp,[],k)';
out = out(:)';

end

function [data_out] = scrambler(data_in,init_state)

if nargin == 1
  init_state = 'CBA63';
end

scrambler_state =de2bi(hex2dec(init_state),23,'left-msb');
data_out = zeros(size(data_in));
for ii = 1:length(data_in)
    data_out(ii) = xor(data_in(ii),xor(scrambler_state(18),scrambler_state(23)));
    scrambler_state(2:end) = scrambler_state(1:end-1);
    scrambler_state(1) = data_out(ii);
end

end

function [data_out] = descrambler(data_in,init_state)

if nargin == 1
  init_state = 'CBA63';
end

scrambler_state =de2bi(hex2dec(init_state),23,'left-msb');
data_out = zeros(size(data_in));
for ii = 1:length(data_in)
    data_out(ii) = xor(data_in(ii),xor(scrambler_state(18),scrambler_state(23)));
    scrambler_state(2:end) = scrambler_state(1:end-1);
    scrambler_state(1) = data_in(ii);
end

end

function [out] = findEVM(a,b)

out = 10*log10(sum((a-b).^2));

end

function [varargout] = spectrumVisualizer(data,fs,nfft,titl,colr,windw)

narginchk(1,6);
nargoutchk(0,2);

data = data(:);
if(nargin == 1)
    fs = 1;
    if(numel(data)<8192)
        nfft    = 2^floor(log2(numel(data)));
    else
        nfft    = 8192;
    end
    titl = 'Power spectral density';
    colr = '-r';
    windw = 'hann';
elseif(nargin == 2)
    if(numel(data)<8192)
        nfft    = 2^floor(log2(numel(data)));
    else
        nfft    = 8192;
    end
    titl = 'Power spectral density';
    colr = '-r';
    windw = 'hann';
elseif(nargin == 3)
    titl = 'Power spectral density';
    colr = '-r';
    windw = 'hann';
elseif(nargin == 4)
    colr = ['-' num2str(colr)];
    windw = 'hann';
elseif(nargin == 5)
    windw = 'hann';
end

if(isempty(fs))
    fs = 1;
end
if(isempty(nfft))
    if(numel(data)<8192)
        nfft    = 2^floor(log2(numel(data)));
    else
        nfft    = 8192;
    end
end
if(isempty(titl))
    titl = 'Power spectral density';
end
if(isempty(colr))
    colr = '-r';
end

if(strcmp(windw,'hann'))
    data = data.*hann(length(data));
elseif(strcmp(windw,'hamming'))
    data = data.*hamming(length(data));
elseif(strcmp(windw,'blackman'))
    data = data.*blackman(length(data));
else
    warning('Window instruction not proper, bypassing the operation');
end

data = data.';

fft_out     = fft(data,nfft)/numel(data);
f_axis      = (-nfft/2+1:nfft/2)*fs/nfft;

abs_fft     = abs(fft_out);
log_abs_fft = 20*log10(abs_fft);
log_abs_fft = fftshift(log_abs_fft);

if nargout == 0
    if(fs ==1)
        x_label = 'Normalized Frequency (-Fs/2 --> Fs/2)';
        title_card = sprintf('%s at Fs = %0.2fHz',titl,fs);
    elseif(fs > 1 && fs <= 1e3)
        x_label = 'Freq(Hz)';
        if(fs == 1e3)
            title_card = sprintf('%s at Fs = %0.2fKHz',titl,fs/1e3);
        else
            title_card = sprintf('%s at Fs = %0.2fHz',titl,fs);
        end
    elseif(fs > 1e3 && fs <= 1e6)
        f_axis  = f_axis/1e3;
        x_label = 'Freq(KHz)';
        if(fs == 1e6)
            title_card = sprintf('%s at Fs = %0.2fMHz',titl,fs/1e6);
        else
            title_card = sprintf('%s at Fs = %0.2fKHz',titl,fs/1e3);
        end
    elseif(fs > 1e6 && fs <= 1e9)
        f_axis  = f_axis/1e6;
        x_label = 'Freq(MHz)';
        if(fs == 1e9)
            title_card = sprintf('%s at Fs = %0.2fGHz',titl,fs/1e9);
        else
            title_card = sprintf('%s at Fs = %0.2fMHz',titl,fs/1e6);
        end
    end

    plot(f_axis,log_abs_fft,colr,'Linewidth',1);
    hold on;grid on;
    x = 0; y = get(gca,'ylim');
    plot([x x],y,'-k','Linewidth',2)

    xlabel(sprintf('%s',x_label));
    ylabel('Magnitude (dB)');
    title(title_card);
elseif nargout == 2
    varargout{1} = f_axis;
    varargout{2} = log_abs_fft;
end
end

function [final_out] = PolyPhaseInterpolator(b,M,inpt)

b    = b(:).';
inpt = inpt(:).';

if(rem(length(b),M)~=0)
    error('Number of Filter taps should be multiple of Interpolation factor');
end
b   = reshape(b,M,[]);

out = zeros(M,length(inpt));
for ii = 1:M
    out(ii,:) = filter(b(ii,:),1,inpt);
end
final_out = out(:);
end

function [final_out] = PolyPhaseDecimator(b,M,inpt)

b    = b(:).';
inpt = inpt(:).';

if(rem(length(b),M)~=0)
    error('Number of Filter taps should be multiple of decimation factor')
end
if(rem(length(inpt),M)~=0)
    error('Number of Input elements should be multiple of decimation factor')
end

out  = zeros(M,length(inpt)/M);
b    = reshape(b,M,[]);
inpt = reshape(inpt,M,[]);
inpt = fliplr(inpt.').'; % flipping the matrix upside down

for ii = 1:M
    out(ii,:) = filter(b(ii,:),1,inpt(ii,:));
end
final_out = sum(out);
end

function out = turboLikeInterleaver(inpt,p,q,j)
% source : 802.22 2011 document 9.6.2

inpt    = inpt(:).';
k       = length(inpt);
pattern = findPattern(k,p,q,j);

out     = inpt(pattern);
end

function pattern = findPattern(k,p,q,j)
% a utility function for turbolike interleaver/deinterleaver

pattern = 1:k;

for itr = 1:j
    for idx = 1:k
        pattern(idx) = mod(k-p+idx+q*p*mod(-idx-p*pattern(idx),k),k);
    end
end

pattern  = pattern + 1;

end

function out = turboLikeDeInterleaver(inpt,p,q,j)
% source : 802.22 2011 document 9.6.2

inpt    = inpt(:).';
k       = length(inpt);
pattern = findPattern(k,p,q,j);

[~,idx] = sort(pattern);
out     = inpt(idx);
end

function symbls = bits2ModSymbols(inpt,M,qam_psk)

if nargin == 1
  M = 4;
  qam_psk = 'psk';
elseif nargin == 2
  qam_psk = 'qam';
end

msgs    = bits2Msg(inpt,log2(M));

if strcmpi(qam_psk,'qam')
  symbls = qammod(msgs,M);
elseif strcmpi(qam_psk,'psk')
  symbls = pskmod(msgs,M);
end

end

function bits = modSymbols2Bits(inpt,M,qam_psk)

if nargin == 1
  M = 4;
  qam_psk = 'psk';
elseif nargin == 2
  qam_psk = 'qam';
end


if strcmpi(qam_psk,'qam')
  msgs = qamdemod(inpt,M);
elseif strcmpi(qam_psk,'psk')
  msgs = pskdemod(inpt,M);
end

bits    = msg2Bits(msgs,log2(M));

end

function encOut = linearBinaryBlockEncoder(msgBits,n,k)

nParityBits = n - k;

if n == 2^nParityBits-1 % hamming code
  encOut  = encode(msgBits,n,k,'hamming/binary');
else
  pol     = cyclpoly(n,k);
  parmat  = cyclgen(n,pol);
  genmat  = gen2par(parmat);
  encOut  = encode(msgBits,n,k,'linear/binary',genmat);
end

end

function decOut = linearBinaryBlockDecoder(msgBits,n,k)

nParityBits = n - k;

if n == 2^nParityBits-1 % hamming code
  decOut  = decode(msgBits,n,k,'hamming/binary');
else
  pol     = cyclpoly(n,k);
  parmat  = cyclgen(n,pol);
  genmat  = gen2par(parmat);
  decOut  = decode(msgBits,n,k,'linear/binary',genmat);
end

end
