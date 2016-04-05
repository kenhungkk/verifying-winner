% Plots the power of the selective test vs. Gupta & Nagel for multinomial
% distribution
% m = number of votes, i.e. X1 + ... + Xn
% n = number of candidates
% alpha = level of the test
% numSample = number of samples to take for simulation purposes
function [] = multinomialPower (m, n, alpha, numSample)
    % delta ranges from 0 to 2
    delta = 0:0.01:2;
    sel = arrayfun(@(x) selectivePower(m, n, x, alpha, numSample), delta);
    gn = gnPower(m, n, delta, alpha, numSample);
    plot(delta, sel, delta, gn, '--');
    legend('Selective', 'Gupta and Nagel', 'Location', 'southeast');
    xlabel('\delta');
    ylabel('Power');
end

% Finds the selective test power for multinomial distribution
% m = number of votes
% n = number of candidates
% distribution is proportional to (exp(delta), 1, ..., 1)
% alpha = level of the test
% numSample = number of samples to take for simulation purposes
function [power] = selectivePower(m, n, delta, alpha, numSample)
    % Probability vector of multinomial distribution
    pi = pr(n, delta);
    % Generates the sample
    X = mnrnd(m, pi, numSample);
    % secX is max(Xj; j>1), i.e. X2
    secX = max(X(:, 2:end), [], 2);
    % nX is X1 + X2
    nX = secX + X(:, 1);
    accepted = nnz(cdf('Binomial', secX, nX, 0.5) * 2 < alpha);
    power = accepted / numSample;
end

% Find the test power of Gutpa and Nagel for multinomial distribution
% m = number of votes
% n = number of candidates
% delta is a vector listing all the deltas we want to plot
% alpha = level of the test
% numSample = number of samples to take for simulation purposes
function [power] = gnPower(m, n, delta, alpha, numSample)
    % Initialization; dStar is the biggest d
    dStar = -Inf;
    for r = 2:n
        % Probability vector (1/r, ..., 1/r)
        pi = ones([1, r]) ./ r;
        % Generates the sample
        X = mnrnd(m, pi, numSample);
        % Takes the maximum of each row
        [maxX, I] = max(X, [], 2);
        % Sets the maximum of each row to -Inf
        X(sub2ind(size(X), 1:numSample, transpose(I))) = -Inf;
        % Takes the maximum of each row again, yielding second largest
        secX = max(X, [], 2);
        % d for a given r is the upper alpha quantile, rounded up
        d = ceil(quantile(maxX - secX, 1 - alpha));
        % Taking maximum over all r
        if dStar < d
            dStar = d;
        end
    end
    % Initialization; power is a vector of powers corresponding to the
    % deltas
    power = zeros(size(delta));
    for i = 1:numel(delta)
        % Probability vector
        pi = pr(n, delta(i));
        % Generates the sample
        X = mnrnd(m, pi, numSample);
        % secX is max(Xj; j>1), i.e. X2
        secX = max(X(:, 2:end), [], 2);
        accepted = nnz(X(:, 1) - dStar > secX);
        power(i) = accepted / numSample;
    end
end

% Returns a probability vector proportional to (exp(delta), 1, ..., 1)
function [pi] = pr(n, delta)
    pi = ones([1, n]);
    pi(1) = exp(delta);
    pi = pi ./ norm(pi, 1);
end