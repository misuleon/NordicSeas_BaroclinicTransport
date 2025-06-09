function [A, amp] = deseason(tim, var, fmt)
% DESEASON Remove seasonal cycle from a time series by subtracting monthly climatology.
%
%   [A, amp] = deseason(tim, var, fmt)
%
%   Inputs:
%       tim - time vector (same length as var)
%             Format depends on `fmt`:
%                1 : MATLAB datenum format
%                2 : Year-decimal format (e.g., 1995.75)
%
%       var - data vector (same size as tim)
%       fmt - time format flag (1 or 2)
%
%   Outputs:
%       A   - deseasoned variable (seasonal cycle removed, long-term mean retained)
%       amp - seasonal amplitude (std. dev. of monthly climatology)
%

    % Get month from input time
    if fmt == 1
        [~, mo, ~] = datevec(tim);
    elseif fmt == 2
        mo = floor(rem(tim, 1) * 12) + 1;
    else
        error('Invalid time format. Use 1 for datenum, 2 for year.decimal');
    end

    % Remove overall mean to isolate seasonal cycle
    valid = ~isnan(var);
    mean_val = mean(var(valid));
    var_anom = var - mean_val;

    % Compute monthly climatology
    monthly_mean = NaN(1, 12);
    for m = 1:12
        monthly_mean(m) = mean(var_anom(mo == m), 'omitnan');
    end

    % Seasonal amplitude (std. dev. of climatology)
    amp = std(monthly_mean, 'omitnan');

    % Subtract seasonal cycle, re-add mean
    A = var - monthly_mean(mo);

    % (Optional debug plot)
    % figure;
    % plot(tim, var, 'k', tim, A, 'r');
    % legend('Original', 'Deseasoned');
end