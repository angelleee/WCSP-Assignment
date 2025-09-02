clear all;
close all;
clc;

%% Parameter Setting

SIR = 0;                      		% signal-to-CCI ratio
IPW = 10^(-SIR/10);                   % interference power
Q = 8; % number of active users
L = 5; % number of paths (fingers)
N = 48; % random code length
trial = 5000;                         % number of Monte Carlo runs
SNR = -20:2:20;                       % signal to noise ratio

%-------------------------------
rMS = zeros(trial,length(SNR));
rZF = zeros(trial,length(SNR));
rRAKE = zeros(trial,length(SNR));

%% Simulation experiments

for jj = 1:trial
    for ii = 1:length(SNR)
        
        NPW = 10^(-SNR(ii)/10);                                       % noise power
        
        %----------------------------fading matrix
        for q = 1:Q
            fad(q,1:L+1) = randn(1,L+1);                              % fading gain
        end
        for qq = 1:Q                                                % fading  normalization
            fad(qq,:) = fad(qq,:)./sqrt(sum(abs(fad(qq,:)).^2));
        end
        
        %----------------------------spreading code
        SP_code = 2*round(rand(Q,N))-1;
        
        %----------------------------perfect channel estimation for H and h1
        H = zeros(N+L,Q);
        for q = 1:Q
            H(:,q) = conv(fad(q,:),SP_code(q,:)).';
        end
        h1 = H(:,1);
        
        %----------------------------correlation matrix generation
        Rin = IPW*H(:,2:Q)*eye(Q-1)*H(:,2:Q)'+NPW*eye(N+L);
        
        %----------------------------SINR for RAKE
        dRAKE = h1;                                               % RAKE receiver
        rRAKE(jj,ii) = norm(dRAKE'*h1)^2/(dRAKE'*Rin*dRAKE);      % RAKE receiver output SINR
        
        %----------------------------SINR for ZF
        dZF = H*inv(H'*H)*[1 zeros(1,Q-1)].';                     % ZF receiver
        rZF(jj,ii) = norm(dZF'*h1)^2/(dZF'*Rin*dZF);          	% ZF receiver output SINR
        
        %----------------------------SINR for MMSE
        Rxx = H*diag([1 IPW*ones(1,Q-1)])*H'+NPW*eye(N+L);
        rxs = h1;
        dMS = inv(Rxx)*rxs;                                       % MMSE detector
        rMS(jj,ii) = norm(dMS'*h1)^2/(dMS'*Rin*dMS);				% MMSE detector output SINR
        
    end
end

rRAKE_AV = sum(rRAKE)/trial;
rZF_AV = sum(rZF)/trial;
rMS_AV = sum(rMS)/trial;

%% Simulation results

H1 = figure;
plot(SNR,10*log10(rMS_AV),'-b+',SNR,10*log10(rZF_AV),'-r+',SNR,10*log10(rRAKE_AV),'-g+');
h = legend('MMSE','ZF','RAKE','location','NorthWest');
set(h,'Fontsize',15);
xlabel('SNRi (dB)');
ylabel('SINRo (dB)');
