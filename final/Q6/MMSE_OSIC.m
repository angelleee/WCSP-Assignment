function [ symbol_MMSE_OSIC ] = MMSE_OSIC(Da_Str,y,H,NPW)
    H_remain = H;
    y_remain = y;
    symbol_MMSE_OSIC = zeros(Da_Str, 1);

    idx_map = 1:Da_Str;

    for i = 1:Da_Str
        
        wmmse_t = (H_remain'*H_remain + size(H_remain,1)*NPW*eye(size(H_remain,2))) \ H_remain';
        G = wmmse_t';

        [~, idx] = min(vecnorm(G, 2, 1));
        wj = G(:,idx);
        zj = wj' * y_remain;
        
        detect = qamdemod(zj*sqrt(2), 4);
        x_hat_rebuild = qammod(detect, 4) / sqrt(2);

        symbol_MMSE_OSIC(idx_map(idx)) = detect;

        y_remain = y_remain - H_remain(:,idx) * x_hat_rebuild;
        H_remain(:,idx) = [];
        idx_map(idx) = [];
    end
end
