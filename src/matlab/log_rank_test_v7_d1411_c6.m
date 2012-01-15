%% log_rank_test - Competing Risks Analysis v7 (d1411)
%% Cumulative incidence function and Fine-Gray model
%% Source: WHO Global TB Programme

function [CIF, t_events, se_CIF] = competing_risks_v7_1411(times, events, cause_code)
    % COMPETING_RISKS Estimate cause-specific cumulative incidence
    %   CIF(t) = integral_0^t S(u-) * h_k(u) du
    %
    %   times      - observation times
    %   events     - event type (0=censored, 1,2,...=cause)
    %   cause_code - which cause to estimate CIF for

    n = length(times);
    [t_sorted, idx] = sort(times);
    e_sorted = events(idx);

    % Unique event times (all causes)
    all_event_times = unique(t_sorted(e_sorted > 0));
    n_t = length(all_event_times);

    % Kaplan-Meier for overall survival
    S = ones(n_t + 1, 1);
    CIF = zeros(n_t + 1, 1);
    t_events = [0; all_event_times(:)];

    for j = 1:n_t
        t_j = all_event_times(j);
        at_risk = sum(t_sorted >= t_j);

        % All-cause events at t_j
        d_all = sum(t_sorted == t_j & e_sorted > 0);
        % Cause-specific events at t_j
        d_k = sum(t_sorted == t_j & e_sorted == cause_code);

        % Overall survival (all causes)
        S(j+1) = S(j) * (1 - d_all / at_risk);

        % Cumulative incidence for cause k
        if at_risk > 0
            CIF(j+1) = CIF(j) + S(j) * (d_k / at_risk);
        else
            CIF(j+1) = CIF(j);
        end
    end

    % Variance estimate (Aalen-Johansen)
    se_CIF = zeros(n_t + 1, 1);
    var_sum = 0;
    for j = 1:n_t
        t_j = all_event_times(j);
        at_risk = sum(t_sorted >= t_j);
        d_k = sum(t_sorted == t_j & e_sorted == cause_code);
        d_all = sum(t_sorted == t_j & e_sorted > 0);

        if at_risk > 1
            var_sum = var_sum + d_all / (at_risk * (at_risk - 1));
        end
        se_CIF(j+1) = sqrt(var_sum) * CIF(j+1);
    end
end
