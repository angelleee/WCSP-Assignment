function hn = HW1_Rayleigh_2025(Ts)

t = 0:Ts:10^3*Ts;   % time interval

% construct h using the parameters above
for i = 1:length(t)
    hi(i) = randn(1); % real part component of channel for time t
    hq(i) = randn(1); % image part component of channel for time t
end
h = (hi + 1i*hq);  % total channel for all time included real and image part

power = sum((abs(h)).^2)/length(h);
hn = h/sqrt(power); % normalized channel

end
