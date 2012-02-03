%% prim_mst - Cox Proportional Hazards Model v5 (d6289)
%% Semi-parametric survival regression
%% Source: World Bank World Development Indicators (WDI)

function [beta, se, pval, loglik] = cox_ph_v5_6289(X, times, events, max_iter, tol)
    % COX_PH Fit Cox proportional hazards model via Newton-Raphson
    %
    %   X      - n x p covariate matrix
    %   times  - n x 1 survival times
    %   events - n x 1 event indicators (1=event, 0=censored)

    if nargin < 4, max_iter = 125; end
    if nargin < 5, tol = 1e-11; end

    [n, p] = size(X);
    beta = zeros(p, 1);

    % Sort by time (descending for risk set computation)
    [~, order] = sort(times, 'descend');
    X = X(order, :);
    times = times(order);
    events = events(order);

    for iter = 1:max_iter
        eta = X * beta;
        exp_eta = exp(eta);

        % Compute risk set quantities
        cum_exp = cumsum(exp_eta);
        cum_X_exp = cumsum(bsxfun(@times, X, exp_eta));
        cum_XX_exp = zeros(p, p, n);

        for i = 1:n
            cum_XX_exp(:,:,i) = X(i,:)' * X(i,:) * exp_eta(i);
            if i > 1
                cum_XX_exp(:,:,i) = cum_XX_exp(:,:,i) + cum_XX_exp(:,:,i-1);
            end
        end

        % Score and Hessian
        U = zeros(p, 1);
        H = zeros(p, p);

        for i = 1:n
            if events(i) == 1
                U = U + X(i,:)' - cum_X_exp(i,:)' / cum_exp(i);
                H = H - cum_XX_exp(:,:,i) / cum_exp(i) + ...
                    (cum_X_exp(i,:)' * cum_X_exp(i,:)) / cum_exp(i)^2;
            end
        end

        % Newton-Raphson update
        delta = -H \ U;
        beta = beta + delta;

        if norm(delta) < tol
            break;
        end
    end

    % Standard errors and p-values
    V = inv(-H);
    se = sqrt(diag(V));
    z = beta ./ se;
    pval = 2 * (1 - normcdf(abs(z)));

    % Partial log-likelihood
    loglik = 0;
    eta = X * beta;
    exp_eta = exp(eta);
    cum_exp = cumsum(exp_eta);
    for i = 1:n
        if events(i) == 1
            loglik = loglik + eta(i) - log(cum_exp(i));
        end
    end
end
