%===================================
% This is the simulator for OFDM.
%===================================
%
%  The function contains a completed SCFDMA transceiver, some parameters need to
%  be defined before running. These parameters are described for instance as the follow
%
%   SP.FFTsize = 128;           % size of the transmitter IFFT and the receiver FFT
%   SP.mod_size = 4;            % specified QPSK transmission
%   SP.inputBlockSize = 128;    % input data block size, i.e. number of sub-carrier within chunk
%   SP.CPsize = 16;             % CP length
%   SP.SNR = [0:5:30];          % simulated SNR range in dB
%   SP.numRun = 10;             % number of simulation iterations
%   SP.channel = pedAchannel;   % defined the channel response
%
%   The final SER will be returned.

function SER = OFDM_system(SP)

SER = zeros(length(SP.SNR),1);
numSymbols = SP.FFTsize;                % number of transmitted OFDM symbols

% initial for M-QPSK symbol generation
BitPerSyml = log2(SP.mod_size);         % calculate bit number per symbol
% SgnlPwr = sum(abs(qammod(0: SP.mod_size-1, SP.mod_size)).^2)/SP.mod_size;
WegVect = pow2([log2(SP.mod_size)-1:-1:0]);

% construct channel
% IID Rayleigh faing channel
tmp_fad = randn(2, length(SP.channel));
raychannel = (tmp_fad(1, :) + 1i*tmp_fad(2, :))/sqrt(2);

% normalize the channel
raychannel = raychannel/sqrt(sum(raychannel.^2));

SP.channel = raychannel.*SP.channel;
H_channel = fft(SP.channel, SP.FFTsize);    % frequency domain version of the channel response


for n = 1:length(SP.SNR)
    
    CurrentSNR = SP.SNR(n);                 % show cuurent status (SNR)
    NPW = 10^(-CurrentSNR/10);      % show noise power status (SNR)
    
    % initialize the error count
    ErrCnt = 0;                             % symbol error count
    
    for k = 1:SP.numRun
        
        correct = 0;
        raw_bits = round(rand(BitPerSyml,SP.inputBlockSize));        % generate N QPSK symbol, each symbol contains n bits
        inputSymbols = qammod(WegVect * raw_bits, SP.mod_size);     % constellation mapping
        inputSamples_ofdm = ifft(inputSymbols, SP.FFTsize);         % convert the signal back to time domain
        
        % add CP
        T_CP_Smples = [inputSamples_ofdm(numSymbols-SP.CPsize+ 1:numSymbols) inputSamples_ofdm];
        
        % propagate through multi-path channel. i.e take the convolution
        R_CP_Smples = conv( T_CP_Smples,SP.channel);
        R_CP_Smples = R_CP_Smples(1:length(T_CP_Smples));
        
        % generate AWGN with appropriate noise power and add it into RxSmples
        Noise = sqrt(NPW/SP.FFTsize)*(randn(1,numSymbols+SP.CPsize) + 1i*randn(1, numSymbols+SP.CPsize));
        R_CP_Smples = R_CP_Smples + Noise;
        
        
        % remove CP
        R_Smples = R_CP_Smples(SP.CPsize+1:numSymbols+SP.CPsize);
        
        % add CFO
        shift = 0:1/SP.FFTsize:(SP.FFTsize-1)/SP.FFTsize;
        R_CFO = exp(1i*2*pi*SP.CFO*shift);
        R_Smples = R_Smples.*R_CFO;
        
        % convert the received signal into frequency domain
        Y_ofdm = fft(R_Smples, SP.FFTsize);
        
        % find the channel response for the subcarriers
        H_eff = H_channel(1:SP.inputBlockSize);
        
        % perform zero forcing equalization in the frequency domain
        X_zf_ofdm = Y_ofdm./H_eff;
        
        % perform hard decision detection
        Det_Symbols = HardDet(X_zf_ofdm);
        
        
        % find and count errors
        for ii = 1:length(Det_Symbols)
            if inputSymbols(ii)-Det_Symbols(ii) == 0
                correct = correct+1;
            end
        end
        ErrCnt = ErrCnt + (length(Det_Symbols)-correct);
        
    end
    
    % calculate the symbol error rate (SER)
    SER(n,:) = ErrCnt / (SP.inputBlockSize*SP.numRun);
    
end


end
