clear all;
close all;
clc;

%% Parameter Setting

Ns = 2000;                          % samples for estimation related parameters
Q=3; % number of active users
L=5; % number of paths (fingers)
N=48; % random code length
trial = 100;                      	% number of Monte Carlo runs

%--------------------------signal power setting
SNR_para = -2:1:4;
Signal_pow = zeros(Q, length(SNR_para));
for q = 1:length(SNR_para)
    Signal_pow(:, q) = (0+5*SNR_para(q)):-6:(-12+5*SNR_para(q));
end

BER = zeros(2, length(SNR_para));

%% Simulation experiments

for SNR_count = 1:length(SNR_para)
    for jj = 1:trial
        
        %--------------------------symbol matrix
        Sa_detect = zeros (Q, Ns);
        Sb_detect = zeros (Q, Ns);
        S = 2*randi([0 1], Q, Ns)-1;
        
        %--------------------------noise matrix
        NPW = 1000;                                       	% noise power
        noise = sqrt(NPW/2)*(randn(N+L, Ns)+1i*randn(N+L, Ns));     % noise
        
        for ii = 1:Ns
            
            %-------------------------fading matrix
            for q = 1:Q
                fad(q, 1:L+1) = randn(1, L+1);                      % fading gain
            end
            for qq = 1:Q                                            % fading normalization
                fad(qq, :) = fad(qq, :)./sqrt(sum(abs(fad(qq, :)).^2));
            end
            
            %----------------------------spreading code
            CODE = 2*randi([0 1], Q, N)-1;
            
            %----------------------------perfect channel estimation for H
            for q = 1:Q
                H(:, q) = conv(fad(q, :),CODE(q, 1:N))';
            end
            
            %----------------------------transmitted symbol
            X1 = H*(NPW*diag([10^(Signal_pow(1, SNR_count)/10), 10^(Signal_pow(2, SNR_count)/10), 10^(Signal_pow(3, SNR_count)/10)])).^0.5*S(:, ii) + noise(:, ii);
            Y1 = X1;
            
            %----------------------------detection SICa
            Z_rake1 = H(:, 1)'*X1;
            if (real(Z_rake1) >= 0)
                Dec_1 = 1;
            else
                Dec_1 = -1;
            end
            
            %----------------------------reconstruction & cancellation
            x_1 = Dec_1*H(:, 1)*(NPW*10^(Signal_pow(1, SNR_count)/10)).^0.5;
            X2 = X1-x_1;
            
            %----------------------------detection
            Z_rake2 = H(:, 2)'*X2;
            if (real(Z_rake2) >= 0)
                Dec_2 = 1;
            else
                Dec_2 = -1;
            end
            
            %----------------------------reconstruction & cancellation
            x_2 = Dec_2*H(:, 2)*(NPW*10^(Signal_pow(2,SNR_count)/10)).^0.5;
            X3 = X2-x_2;
            
            %----------------------------detection
            Z_rake3 = H(:, 3)'*X3;
            if (real(Z_rake3) >= 0)
                Dec_3 = 1;
            else
                Dec_3 = -1;
            end
            
            %----------------------------detection
            Sa_detect(:, ii) = [Dec_1 Dec_2 Dec_3];
            
            %----------------------------detection SICb
            Z_rake1 = H(:, 3)'*Y1;
            if (real(Z_rake1) >= 0)
                Dec_1 = 1;
            else
                Dec_1 = -1;
            end
            
            %----------------------------reconstruction & cancellation
            y_1 = Dec_1*H(:, 3)*(NPW*10^(Signal_pow(3, SNR_count)/10)).^0.5;
            Y2 = Y1-y_1;
            
            %----------------------------detection
            Z_rake2 = H(:, 2)'*Y2;
            if (real(Z_rake2) >= 0)
                Dec_2 = 1;
            else
                Dec_2 = -1;
            end
            
            %----------------------------reconstruction & cancellation
            y_2 = Dec_2*H(:, 2)*(NPW*10^(Signal_pow(2, SNR_count)/10)).^0.5;
            Y3 = Y2-y_2;
            
            %----------------------------detection
            Z_rake3 = H(:, 1)'*Y3;
            if (real(Z_rake3) >= 0)
                Dec_3 = 1;
            else
                Dec_3 = -1;
            end
            
            %----------------------------detection
            Sb_detect(:, ii) = [Dec_3 Dec_2 Dec_1];
            
        end
        
        %----------------------------BER calculation
        error = 0;
        for k = 1:Ns
            for n = 1:Q
                if (S(n, k) ~= Sa_detect(n, k))
                    error = error+1;
                end
            end
        end
        BER(1, SNR_count) = error/(Q*Ns)+BER(1, SNR_count);
        error = 0;
        for k = 1:Ns
            for n = 1:Q
                if (S(n, k) ~= Sb_detect(n, k))
                    error = error+1;
                end
            end
        end
        BER(2, SNR_count) = error/(Q*Ns)+BER(2, SNR_count);
        
    end
end

BER = BER/trial;

%% Simulation results

H = figure;
semilogy(Signal_pow(1, :), BER(1, :), '-o', Signal_pow(1, :), BER(2, :), '-x');
% plot(Signal_pow(1, :), BER(1, :), '-o', Signal_pow(1, :), BER(2, :), '-x');
grid on;
legend('SIC_a', 'SIC_b');
xlabel('SNR (dB)');
ylabel('BER');
