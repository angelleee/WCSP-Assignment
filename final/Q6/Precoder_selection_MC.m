function [ precoder ] = Precoder_selection_MC(codebook,H,Da_Str,noise_power)
    precoder_num = size(codebook, 3);
    max_capacity = -inf;
    precoder = [];

    for i = 1 : precoder_num
        F = codebook(:,:,i);
        Heff = H * F;
        capacity = log2(det( eye(Da_Str) + (1/noise_power) * (Heff' * Heff) ));
        if capacity > max_capacity
            max_capacity = capacity;
            precoder = F;
        end
    end
end
