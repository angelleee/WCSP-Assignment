function [ symbol_K_best ] = K_best(Da_Str,y,H,Q,R) 
    K = 6;
    const = qammod(0:3, 4) / sqrt(2);
    y_hat = Q' * y;
    paths = struct('symbols', [], 'metric', 0);
    
    for level = Da_Str:-1:1
        % build the whole map per level
        new_paths = [];
        for i = 1:length(paths)
            path = paths(i);
            for s = const
                sym = [s; path.symbols];
                r = R(level, level:Da_Str) * sym;
                metric = abs(y_hat(level) - r(1))^2 + path.metric;
                new_paths(end+1).symbols = sym;
                new_paths(end).metric = metric;
            end
        end

        % keep K nodes
        [~, idx] = sort([new_paths.metric]);
        paths = new_paths(idx(1:min(K, end)));
    end
    
    detected_symbols = paths(1).symbols;
    symbol_K_best = qamdemod(detected_symbols * sqrt(2), 4);
end
