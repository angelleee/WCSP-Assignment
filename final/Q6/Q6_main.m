clc
close all
clear
tic
Nt = 4;
Nr = 4;

Da_Str = 4; % number of data streams
signal_power = 1;
SNR = -10:5:20;
noise_power = signal_power./(10.^(SNR./10));
channel_rho = 1; % channel correlation parameter % 0/0.5/1
precoding_method = 2; % 1: codebook-based; 2: SVD 
channel_realization = 100;
trial = 1000;
codebook = Codebook_Generation;  % 4 x 4 x 8

Pe_MMSE = zeros(length(SNR), channel_realization, trial);
Pe_MMSE_OSIC = zeros(length(SNR), channel_realization, trial);
Pe_Kbest = zeros(length(SNR), channel_realization, trial);
Pe_ML = zeros(length(SNR), channel_realization, trial);

for nn = 1:length(SNR)  % you can use parfor here to accelerate your code
    for hh = 1:channel_realization  
        H = Channel_Construction(channel_rho,Nt,Nr);
        %% Symbol Generation
        symbol_d = randi([0, 3], [Da_Str, 1]);  % QPSK: 0~3
        data_symbol = qammod(symbol_d, 4) / sqrt(2);
        %% Precoder Selection
        switch(precoding_method)
            case 1  % codebook-based precoding
                Fopt = Precoder_selection_MC(codebook,H,Da_Str,noise_power(nn));
            case 2  % SVD precoding
                [U, S, V] = svd(H);
                Fopt = V(:, 1:Da_Str);
        end
        [Q,R] = qr(H*Fopt);
        
        for tt = 1:trial
            n = sqrt(noise_power(nn)) * randn(Nr,1, 'like', 1j);
            %% Received Signal Generation
            y = H * Fopt * data_symbol + n;
            
            %% MMSE, MMSE-SIC, K-best SD,ML
            symbol_MMSE = MMSE(Da_Str,y,H*Fopt,noise_power(nn));
            symbol_MMSE_OSIC = MMSE_OSIC(Da_Str,y,H*Fopt,noise_power(nn));
            symbol_K_best = K_best(Da_Str,y,H*Fopt,Q,R);
            symbol_ML = ML(Da_Str,y,H*Fopt,Q,R,noise_power(nn));
            
            %% Symbol DetectionPe_Kbest
            [~, error_MMSE] = symerr(symbol_d, symbol_MMSE);
            Pe_MMSE(nn, hh, tt) = error_MMSE;
            
            [~, error_MMSE_OSIC] = symerr(symbol_d, symbol_MMSE_OSIC);
            Pe_MMSE_OSIC(nn, hh, tt) = error_MMSE_OSIC;
            
            [~, error_Kbest] = symerr(symbol_d, symbol_K_best);
            Pe_Kbest(nn, hh, tt) = error_Kbest;

            [~, error_ML] = symerr(symbol_d, symbol_ML);
            Pe_ML(nn, hh, tt) = error_ML;
            
        end
    end
end
%% Calculate Average SER
average_SER_MMSE = mean(Pe_MMSE, [2 3]);
average_SER_MMSE_OSIC = mean(Pe_MMSE_OSIC, [2, 3]);
average_SER_K_best = mean(Pe_Kbest, [2, 3]);
average_SER_ML = mean(Pe_ML, [2, 3]);
%%
semilogy(SNR,average_SER_MMSE,'ro-',SNR,average_SER_MMSE_OSIC,'ks-',SNR,average_SER_K_best,'bv-',SNR,average_SER_ML,'kd-')
xlabel('SNR')
ylabel('average SER')
ylim([10^(-5) 10^0]);
legend('MMSE','MMSE-OSIC','K-best SD','ML','location','SouthWest')
title(sprintf('(a) £l = %g', channel_rho)) % a,b,c,d
grid on
toc