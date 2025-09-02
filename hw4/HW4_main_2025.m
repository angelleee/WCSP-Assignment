% Comparison of ABBA and D-STTD
clc
clear
close all
% ======= parameter setting =========
% range of SNR
SNR = -15:3:15;

Pe_abba = zeros(1, length(SNR));        % define BER for ABBA encoding
Pe_dsttd = zeros(1, length(SNR));       % define BER for dsttd encoding

mod_size = 4;
BitPerSyml = log2(mod_size);
WegVect = pow2([log2(mod_size)-1:-1:0]);

trial_loop = [2000 2000 2000 2000 2000 2000 5000 5000 10000 20000 40000];

% ====== Monte Carlo process ===============
for snr = 1:length(SNR)
        Trial = trial_loop(snr);
    for t = 1:Trial
        
        % ----- generate data -------
        N = 4;  % 4 symbols are transmitted over 4 symbol times for ABBA
                % 4 symboles are transmitted over 2 symbol times for D-STTD
        raw_bits = randi([0 1], BitPerSyml, N);
        data = WegVect * raw_bits;
        inputSymbols = qammod(data, mod_size);

        % -------- transmit data symbol------------------
        tx_MQAM = 1/sqrt(2)*inputSymbols.';
        
        % --------- channel model & received signal model -----
        NPW = 10^(-SNR(snr)/10);    % noise power
        noise = sqrt(NPW/2)*(randn(N, 1)+1i*randn(N, 1));
        
        % equivalent channel model for ABBA
       
        h_abba = 1/((2*4)^0.5)*(randn(1,4)+1i*randn(1,4));
        H_abba = [h_abba(1),h_abba(2),h_abba(3),h_abba(4);...
            conj(h_abba(2)),-conj(h_abba(1)),conj(h_abba(4)),-conj(h_abba(3));...
            h_abba(3),h_abba(4),h_abba(1),h_abba(2);...
            conj(h_abba(4)),-conj(h_abba(3)),conj(h_abba(2)),-conj(h_abba(1))];
        rx_sig_abba = H_abba*tx_MQAM + noise;
        
        % equivalent channel model for D-STTD
        
        h1 = 1/((2*4)^0.5)*(randn(1,4)+1i*randn(1,4));
        h2 = 1/((2*4)^0.5)*(randn(1,4)+1i*randn(1,4));
        H1=[h1;...
            conj(h1(2)),-conj(h1(1)),conj(h1(4)),-conj(h1(3))];
        H2=[h2;...
            conj(h2(2)),-conj(h2(1)),conj(h2(4)),-conj(h2(3))];
        H_dsttd=[H1;H2];
        rx_sig_dsttd = H_dsttd*tx_MQAM + noise;
        
        
        % ======= decoding process =============
        % ---- ABBA -------
         data_abba = MMSE_OSIC_2025(H_abba, rx_sig_abba, N, NPW, mod_size, 2);
         
         % calculate BER
         [~, error_abba] = biterr(data, data_abba);
         Pe_abba(snr) = Pe_abba(snr) + error_abba;
         
      	% ---- D-STTD -------
         data_dsttd = MMSE_OSIC_2025(H_dsttd, rx_sig_dsttd, N, NPW, mod_size, 2);
         
         % calculate BER
         [~, error_dsttd] = biterr(data, data_dsttd);
         Pe_dsttd(snr) = Pe_dsttd(snr) + error_dsttd;
        
   
    end % end of Trial
    Pe_abba(snr) = Pe_abba(snr)/Trial;
    Pe_dsttd(snr) = Pe_dsttd(snr)/Trial;
end  % end of snr
 
semilogy(SNR,Pe_abba,'b-x')
hold on
semilogy(SNR,Pe_dsttd,'r-o')
legend('ABBA (4x1)','D-STTD (4x2)')
xlabel('SNR (dB)')
ylabel('BER')
