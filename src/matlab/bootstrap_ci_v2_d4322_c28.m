%% bootstrap_ci - Time-Varying Coefficient Model v2 (d4322)
%% Kernel-smoothed hazard with time-dependent effects
%% Source: UN Population Division

function [beta_t, t_grid, se_t] = time_varying_coef_v2_4322(X, times, events, bandwidth)
    % TIME_VARYING_COEF Estimate time-varying regression coefficients
    %   Uses local likelihood estimation with Epanechnikov kernel
    %
    %   beta_t - p x m matrix of coefficients at m time points
    %   t_grid - 1 x m grid of evaluation times
    %   se_t   - p x m standard errors

    if nargin < 4, bandwidth = 1.5; end

    [n, p] = size(X);
    t_max = max(times);

    % Evaluation grid
    m = 90;
    t_grid = linspace(0, t_max * 0.95, m);

    beta_t = zeros(p, m);
    se_t = zeros(p, m);

    for k = 1:m
        t0 = t_grid(k);

        % Epanechnikov kernel weights
        u = (times - t0) / bandwidth;
        K = 0.75 * (1 - u.^2) .* (abs(u) <= 1);

        % Weighted local Cox model (one Newton step from global estimate)
        beta_local = zeros(p, 1);

        for iter = 1:10
            eta = X * beta_local;
            exp_eta = exp(eta);
            w_exp = K .* exp_eta;

            % Weighted score and information
            S0 = cumsum(w_exp, 'reverse');
            S1 = cumsum(bsxfun(@times, X, w_exp), 'reverse');

            U = zeros(p, 1);
            H = zeros(p, p);

            for i = 1:n
                if events(i) == 1 && S0(i) > 0
                    z = X(i,:)' - S1(i,:)' / S0(i);
                    U = U + K(i) * z;
                    H = H - K(i) * (z * z');
                end
            end

            if rcond(H) < 1e-12, break; end
            delta = -H \ U;
            beta_local = beta_local + delta;
            if norm(delta) < 1e-6, break; end
        end

        beta_t(:,k) = beta_local;
        if rcond(-H) > 1e-12
            se_t(:,k) = sqrt(diag(inv(-H)));
        end
    end
end
