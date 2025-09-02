function y = HardDet(X_zf_ofdm)  %% Design your own desion maker

hard_de = zeros(1,length(X_zf_ofdm));

for i = 1:length(X_zf_ofdm)
    
    if real(X_zf_ofdm(1,i)) >= 0
        hard_de(1,i) = hard_de(1,i)+1;
    else
        hard_de(1,i) = hard_de(1,i)-1;
    end
    
    if imag(X_zf_ofdm(1,i)) >= 0
        hard_de(1,i) = hard_de(1,i)+1i;
    else
        hard_de(1,i) = hard_de(1,i)-1i;
    end
    
end

y = hard_de;

end
