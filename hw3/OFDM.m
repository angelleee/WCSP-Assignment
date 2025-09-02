clear all;
close all;

CFOsize=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0];

SP.FFTsize = 256;           % size of the transmitter IFFT and the receiver FFT
SP.mod_size = 4;            % adopt QPSK modulation
SP.inputBlockSize = SP.FFTsize;    % input data block size
SP.CPsize = 16;             % CP length
SP.SNR = [0:5:30];            % simulated SNR range in dB
SP.numRun = 10;             % number of simulation iterations
%SP.CFO = 0;                 % CFO is zero
SP.channel = ones(1,10);     % initial channel
Packet_numRun = 1000;       % number of packets

SER_ = zeros(length(SP.SNR), length(CFOsize));

for i = 1:length(CFOsize)
    SP.CFO = CFOsize(i);             % CP length
    
    SER = 0;
    % run the simulation for OFDM with different channel
    for n = 1:1:Packet_numRun
        SER = OFDM_system(SP) + SER;
    end
    SER = SER / Packet_numRun;

    SER_(:,i) = SER;
end



H1=figure(1);
semilogy(SP.SNR, SER_(:,1), SP.SNR, SER_(:,2), SP.SNR, SER_(:,3), SP.SNR, SER_(:,4), SP.SNR, SER_(:,5), SP.SNR, SER_(:,6), SP.SNR, SER_(:,7));
legend('CFO = 0','CFO = 0.1','CFO = 0.2','CFO = 0.3','CFO = 0.4','CFO = 0.5','CFO = 1');
xlabel('SNR (dB)');
ylabel('SER');
grid on;
title('QPSK OFDM with ZF equalizer under different CFO');
%saveas(H1,'OFDM_normal.jpg')
%save OFDM                   % Save the simulation results
