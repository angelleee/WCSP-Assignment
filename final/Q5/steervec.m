function a = steervec(N, theta)

    a = exp(1j * pi * (0 : N-1).' * sind(theta));

end