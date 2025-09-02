% clear all;
close all;
clc;

%% Parameter setting
Ns = 2000; % samples for estimation related parameters
Q = 3;       % number of active users
L = 5;       % number of paths (fingers)
N = 48;       % random code length
trial = 100; % number of Monte Carlo runs

SNR = -10:5:20; % signal-to-noise ratio

BER = zeros(length(SNR), 1);

Amp = diag(sqrt(10.^([0 -6 -12] / 10)));

%% Simulation experiments

for kk = 1:length(SNR)

    for jj = 1:trial

        %-----------------symbol matrix----------------------
        S_detect = zeros (Q, Ns);
        S = 2 * randi([0 1], Q, Ns) - 1;

        %-----------------noise matrix----------------------
        NPW = 10^(-SNR(kk) / 10); % noise power
        noise = sqrt(NPW / 2) * (randn(N + L, Ns) + 1i * randn(N + L, Ns)); % noise

        for ii = 1:Ns

            %------------fading matrix----------------------
            fad = randn(Q, L + 1); % fading gain
            fad = fad ./ sqrt(sum(abs(fad).^2, 2)); % fading normalization

            %------------------spreading code------------------------
            SP_code = 2 * round(rand(Q, N)) - 1; % spreading code

            %------------------perfect channel estimation for H------
            H = zeros(N + L, Q);

            for q = 1:Q
                H(:, q) = conv(fad(q, :), SP_code(q, :)).';
            end

            %------------------transmitted symbol--------------------
            x = H * Amp * S(:, ii) + noise(:, ii);

            %-------------------detection PIC------------------------
            y = x;
            x_regen = zeros(N + L, Q);
            S_tmp_detect = zeros(Q, 1);

            for q = 1:Q
                % RAKE
                d_RAKE = H(:, q);
                z = d_RAKE' * y;
                % decision
                S_tmp_detect(q) = (real(z) >= 0) * 2 - 1;
                % regeneration
                x_regen(:, q) = H(:, q) * Amp(q, q) * S_tmp_detect(q);
            end

            for q = 1:Q
                % substraction
                x_regen_sum = sum(x_regen(:, [1:(q - 1) (q + 1):end]), 2);
                y_sub = y - x_regen_sum;
                % RAKE
                d_RAKE = H(:, q);
                z = d_RAKE' * y_sub;
                % decision
                S_detect(q, ii) = (real(z) >= 0) * 2 - 1;
            end

        end

        %----------------------------BER calculation----------------------
        error = sum(S ~= S_detect, 'all');
        BER(kk) = error / (Q * Ns) + BER(kk);

    end

end

BER = BER / trial;

%% Simulation results

H = figure;

semilogy(SNR, BER, '-bo');
grid on;
legend('PIC');
xlabel('SNR (dB)');
ylabel('BER');
