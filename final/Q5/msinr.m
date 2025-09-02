function w = msinr(x, a0, sp)
    in = x - a0 * sp;
    Rin = (in * in') / size(in, 2);
    w = (Rin \ a0) / (a0' * (Rin \ a0));
end