function [ symbol_ML ] = ML(Da_Str,y,H,Q,R,NPW)
    const = qammod(0:3, 4) / sqrt(2);
    X_all = combvec(const, const, const, const);
    y_hat = Q' * y;
    min_metric = inf;
    best_x = zeros(Da_Str, 1);

    for i = 1:size(X_all, 2)
        x_try = X_all(:,i);
        metric = norm(y_hat - R * x_try)^2;

        if metric < min_metric
            min_metric = metric;
            best_x = x_try;
        end
    end
    
    symbol_ML = qamdemod(best_x * sqrt(2), 4);
end