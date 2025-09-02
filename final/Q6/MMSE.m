function [ symbol_MMSE ] = MMSE(Da_Str,y,H,NPW)
    wmmse_t = (H'*H + size(H,1)*NPW*eye(Da_Str)) \ H';
    x_hat = wmmse_t * y;
    symbol_MMSE = qamdemod(x_hat*sqrt(2), 4);
end