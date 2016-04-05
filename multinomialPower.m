function [] = multinomialPower (m, n, alpha, numSample)
    delta = 0:0.01:2;
    gn = gnPower(m, n, delta, alpha, numSample);
    sel = arrayfun(@(x) selectivePower(m, n, x, alpha, numSample), delta);
    plot(delta, sel, delta, gn, '--');
    legend('Selective', 'Gupta and Nagel', 'Location', 'southeast');
    xlabel('\delta');
    ylabel('Power');
end

function [power] = selectivePower(m, n, delta, alpha, numSample)
    pi = pr(n, delta);
    X = mnrnd(m, pi, numSample);
    secX = max(X(:, 2:end), [], 2);
    nX = secX + X(:, 1);
    accepted = nnz(cdf('Binomial', secX, nX, 0.5) * 2 < alpha);
    power = accepted / numSample;
end

function [power] = gnPower(m, n, delta, alpha, numSample)
    dStar = -Inf;
    for r = 2:n
        pi = ones([1, r]) ./ r;
        X = mnrnd(m, pi, numSample);
        [maxX, I] = max(X, [], 2);
        X(sub2ind(size(X), 1:numSample, transpose(I))) = -Inf;
        secX = max(X, [], 2);
        diffX = sort(maxX - secX);
        d = ceil(quantile(diffX, 1 - alpha));
        if dStar < d
            dStar = d;
        end
    end
    power = zeros(size(delta));
    for i = 1:numel(delta)
        pi = pr(n, delta(i));
        X = mnrnd(m, pi, numSample);
        secX = max(X(:, 2:end), [], 2);
        accepted = nnz(X(:, 1) - dStar > secX);
        power(i) = accepted / numSample;
    end
end

function [pi] = pr(n, delta)
    pi = ones([1, n]);
    pi(1) = exp(delta);
    pi = pi ./ norm(pi, 1);
end