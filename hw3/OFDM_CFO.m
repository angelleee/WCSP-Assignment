clear all;
close all;

SP.FFTsize = 256;           % size of the transmitter IFFT and the receiver FFT
SP.mod_size = 4;            % adopt QPSK modulation
SP.inputBlockSize = SP.FFTsize;    % input data block size
SP.CPsize = 16;             % CP length.
SP.SNR = [0:5:30];            % simulated SNR range in dB
SP.numRun = 10;             % number of simulation iterations
SP.channel = ones(1,10);     % initial channel
Packet_numRun = 1000;       % number of packets

% set CFO
CFO = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 1];
SER = zeros(length(SP.SNR), length(CFO));

% run the simulation for OFDM with different channel
for k = 1:length(CFO)
    SP.CFO = CFO(k);
    for n = 1:Packet_numRun
        SER(:, k) = OFDM_system(SP) + SER(:, k);
    end
    SER(:, k) = SER(:, k) / Packet_numRun;
end

H1=figure(1);
semilogy( SP.SNR, SER(:, 7), '-v',SP.SNR, SER(:, 6), '-*', SP.SNR, SER(:, 5), '-s', SP.SNR, SER(:, 4), '-^', SP.SNR, SER(:, 3) ,'-x', SP.SNR, SER(:, 2), '-+',SP.SNR, SER(:, 1), '-o' );
legend('CFO = 1','CFO = 0.5','CFO = 0.4','CFO = 0.3','CFO = 0.2','CFO = 0.1','CFO = 0','location','SouthWest');
xlabel('SNR(dB)');
ylabel('SER');
grid on;
title('QPSK OFDM with ZF equalizer');
% saveas(H1,'OFDM_CFO.jpg');
