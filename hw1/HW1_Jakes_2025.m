% hn = Jakes(M,f_max,Ts)
% It's a subroutine which is about to generate
% the Doppler spread fading channel based on Jakes model
%
% Input:
%       M: a parameter ralated to no. of paths
%          (no. of path N = 4*M + 2)
%       f_max: maximum Doppler shift (Hz)
%       Ts: sampling rate (sec)
% Output:
%       hn = normalized channel in the interval, [0,10^3*Ts]

function hn = HW1_Jakes_2025(M,f_max,Ts)

N = 4*M + 2;    % no. of paths
alpha = 0;      % a parameter ralated to the inner product properties of the output channel
n = 1:1:M;
t = 0:Ts:10^3*Ts;   % time interval
beta = n*pi/M;

% construct h using the parameters above
fn = f_max*cos(2*pi*n/N); % oscillation frequency
for i = 1:length(t)
    hi(i) = 2*cos(beta)*cos(2*pi*fn*t(i))'+sqrt(2)*cos(alpha)*cos(2*pi*f_max*t(i)); % real part component of channel for time t
    hq(i) = 2*sin(beta)*cos(2*pi*fn*t(i))'+sqrt(2)*sin(alpha)*cos(2*pi*f_max*t(i)); % image part component of channel for time t
end
h = (hi + 1i*hq);  % total channel for all time included real and image part

power = sum((abs(h)).^2)/length(h);
hn = h/sqrt(power); % normalized channel

end
