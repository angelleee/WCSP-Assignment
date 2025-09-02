function w = mmse(x, sp)
    Rxx = (x * x') / size(x, 2);
    rxs = (x * sp') / size(x, 2);
    w = Rxx \ rxs;
end