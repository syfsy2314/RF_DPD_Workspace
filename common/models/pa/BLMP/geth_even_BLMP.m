function H2 = geth_even_BLMP(x, K, Q, h)
%GETH_EVEN_BLMP Construct BLMP basis functions and apply the feedback FIR.

    validateattributes(x, {'numeric'}, {'vector', 'nonempty', 'finite'}, ...
        mfilename, 'x');
    validateattributes(h, {'numeric'}, {'vector', 'nonempty', 'finite'}, ...
        mfilename, 'h');

    fullBandBasis = geth_even_1d(x(:), K, Q);
    H2 = filter(h(:), 1, fullBandBasis);

% Original direct implementation retained for reference.
%{
    s = K * (Q + 1);
    H2 = zeros(length(x), s);
    P = length(h);
    t = 1;

    for k = 1:K
        for q = 0:Q
            for p = 0:(P - 1)
                x1 = [zeros(q + p, 1); x(1:(end - q - p))];
                H2(:, K * q + k) = H2(:, K * q + k) + ...
                    (abs(x1).^(k - 1)) .* x1 .* h(p + 1);
                % H2(:, t) = H2(:, t) + ...
                %     (abs(x1).^(k - 1)) .* x1 .* h(p + 1);
            end
            t = t + 1;
        end
    end

    % Alternative two-stage implementation from the original file:
    % H = zeros(length(x), s);
    % n = 1;
    % for q = 0:Q
    %     x1 = [zeros(q, 1); x(1:(end - q))];
    %     for k = 1:K
    %         H(:, n) = H(:, n) + (abs(x1).^(k - 1)) .* x1;
    %         n = n + 1;
    %     end
    % end
    % filterOrder = length(h) - 1;
    % HH = zeros(size(H));
    % for nn = 0:filterOrder
    %     H_temp = [zeros(nn, s); H(1:(end - nn), :)];
    %     HH = h(nn + 1) .* H_temp + HH;
    % end
%}
end
