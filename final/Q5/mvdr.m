function w = mvdr(x, a0)
    Rxx = (x * x') / size(x, 2);
    w = (Rxx \ a0) / (a0' * (Rxx \ a0));
end