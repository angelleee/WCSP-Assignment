function [ H ] = Channel_Construction(rho,Nt,Nr)

% rho = 0: uncorrelated channel
% rho = 0.5: medium-correlation channel
% rho = 1: fully correlated channel

Hw = (1/sqrt(2))*complex(randn(Nr,Nt),randn(Nr,Nt));

Rt = zeros(Nt);
for ii = 1:Nt
    for jj = 1:Nt
        Rt(ii,jj) = rho^(abs(ii-jj));
    end
end

Rr = zeros(Nr);
for ii = 1:Nr
    for jj = 1:Nr
        Rr(ii,jj) = rho^(abs(ii-jj));
    end
end

H = Rr^(1/2)*Hw*Rt^(1/2);
H = H/sqrt(sum(sum(abs(H).^2)))*sqrt(Nt*Nr);

end
