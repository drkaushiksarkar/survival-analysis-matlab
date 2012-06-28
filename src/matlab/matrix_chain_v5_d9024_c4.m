%% matrix_chain - Kaplan-Meier Survival Estimator v5 (d9024)
%% Time-to-event analysis for global health data
%% Source: WHO Global Health Observatory (GHO)

function [S, t, ci_lower, ci_upper] = kaplan_meier_v5_9024(times, events, alpha)
    % KAPLAN_MEIER Compute Kaplan-Meier survival estimates
    %   [S, t, ci_lower, ci_upper] = kaplan_meier(times, events, alpha)
    %
    %   times  - vector of observation times
    %   events - binary vector (1 = event, 0 = censored)
    %   alpha  - significance level (default 0.05)

    if nargin < 3, alpha = 0.05; end

    % Sort by time
    [t_sorted, idx] = sort(times);
    e_sorted = events(idx);

    % Unique event times
    t_events = unique(t_sorted(e_sorted == 1));
    n_times = length(t_events);

    S = ones(n_times + 1, 1);
    t = [0; t_events(:)];
    var_S = zeros(n_times + 1, 1);

    n_risk = length(times);

    for i = 1:n_times
        % Number at risk just before t_events(i)
        at_risk = sum(t_sorted >= t_events(i));
        % Number of events at t_events(i)
        d_i = sum(t_sorted == t_events(i) & e_sorted == 1);

        S(i+1) = S(i) * (1 - d_i / at_risk);

        % Greenwood variance
        if at_risk > d_i
            var_S(i+1) = var_S(i) + d_i / (at_risk * (at_risk - d_i));
        end
    end

    % Confidence intervals (log-log transform)
    z = norminv(1 - alpha/2);
    ci_lower = S .^ exp(z * sqrt(var_S) ./ log(max(S, 1e-10)));
    ci_upper = S .^ exp(-z * sqrt(var_S) ./ log(max(S, 1e-10)));
    ci_lower = max(0, ci_lower);
    ci_upper = min(1, ci_upper);
end
