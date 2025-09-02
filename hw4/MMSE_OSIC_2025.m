function [decide_data] = MMSE_OSIC_2025(H, in, P, NPW, M, L)
% H : is the channel matrix
% in : is the input data
% P : is the number of data streams
% NPW : is the noise power
% 4QAM used: M=4, L=2
F = H'*H;
in_match = H'*in;
decide_data = zeros(1,P);
dec_table = 1:P; % a index table helps make decision in tx_order

for k = 1:P
    G = inv(F*F' + NPW*F)*F;
    [~, index] = sort(diag(G'*G));
    
    Z=G(:,index(1))'*in_match;
    decide = qamdemod(Z,M);
    decide_data(dec_table(index(1))) = decide;
    dec = qammod(decide,M);
    in=in-H(:,index(1))*1/(L^0.5)*dec;
    
    dec_table(index(1)) = [];
    H(:,index(1)) = [];
    in_match = H'*in;
    F = H'*H;
end

end