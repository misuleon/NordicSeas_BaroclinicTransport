function [yrs, Trnsprt, Trns_uncert] = compute_baroclinic_transport_NS_RT()
% COMPUTE_BAROCLINIC_TRANSPORT_NS_RT
%   Estimate annual baroclinic transport into the Nordic Seas
%   hydrographic profiles in the Norwegian Sea and Rockall Trough.
%   Returns transport time series and associated uncertainties.
%
%   Output:
%       yrs         - vector of years [1900–2023]
%       Trnsprt     - baroclinic transport anomaly (Sv)
%       Trns_uncert - uncertainty in transport (Sv)  

    %% Load data
    load South_Nor_Sea_to_2022_500
    xb = -[2 5 6.5 6.5 2 2]; yb = [63.55 63.55 64 66 66 63.55];
    J = inside3(lon, lat, xb, yb); J = ~isnan(J);
    date_NS = stn_date(J);
    temp_NS = temp(:,J); dyn_NS = dyn(:,J);
    [date_NS, I] = sort(date_NS);
    temp_NS = temp_NS(:,I); dyn_NS = dyn_NS(:,I);
    yrNS = year(date_NS);

    load Rockall_to_2022_500
    xb = -[10.5 14.5 13.3 12.5 12.5 12.8 13.0 10 10 9.5 9.5 10.5];
    yb = [55.5 55.5 57.0 57.4 57.8 58 58.5 58.5 57.5 57 56 55];
    J = inside3(lon, lat, xb, yb); J = ~isnan(J);
    date_RT = stn_date(J);
    temp_RT = temp(:,J); dyn_RT = dyn(:,J);
    [date_RT, I] = sort(date_RT);
    temp_RT = temp_RT(:,I); dyn_RT = dyn_RT(:,I);
    yrRT = year(date_RT);

    %% Deseason dynamic height
    ti_NS = datenum(date_NS) - datenum(1950,1,1); ti_NS = ti_NS / 365 + 1950;
    ti_RT = datenum(date_RT) - datenum(1950,1,1); ti_RT = ti_RT / 365 + 1950;
    for i = 1:length(dpth)
        dyn_NS(i,:) = deseason(ti_NS, dyn_NS(i,:), 2);
        dyn_RT(i,:) = deseason(ti_RT, dyn_RT(i,:), 2);
    end

    %% Reference to 500 m (index 25)
    dyn_NS_ref = dyn_NS(25,:) - dyn_NS;
    dyn_RT_ref = dyn_RT(25,:) - dyn_RT;
    PE_NS = sum_Z(dpth, dyn_NS_ref);
    PE_RT = sum_Z(dpth, dyn_RT_ref);

    % Remove Rockall 1949 outliers
    yrRT([17,18]) = [];
    PE_RT([17,18]) = [];

    %% Annual means and uncertainties
    yrSt = 1900; yrEn = 2023; N_yrs = yrEn - yrSt + 1;
    for i = 1:N_yrs
        y = yrSt + i - 1;
        k = yrNS == y;
        mean_PE_NS(i) = mean(PE_NS(k), 'omitnan');
        std_PE_NS(i)  = std(PE_NS(k), 'omitnan');
        N_NS(i)       = sum(k);

        j = yrRT == y;
        mean_PE_RT(i) = mean(PE_RT(j), 'omitnan');
        std_PE_RT(i)  = std(PE_RT(j), 'omitnan');
        N_RT(i)       = sum(j);
    end

    % Align 1936 NS to RT
    mean_PE_NS(39) = mean_PE_NS(37);
    std_PE_NS(39)  = std_PE_NS(37);

    %% Transport and uncertainty
    for i = 1:N_yrs
        Trnsprt(i) = mean_PE_RT(i) - mean_PE_NS(i);
        Trns_uncert(i) = sqrt(std_PE_NS(i)^2 + std_PE_RT(i)^2) / sqrt(N_NS(i) + N_RT(i));
    end

    eff = 4 * pi / 86400 * sind(60);
    Trnsprt = Trnsprt * 10 / eff / 1e6;
    Trns_uncert = Trns_uncert * 10 / eff / 1e6;
    yrs = yrSt:yrEn;

    %% Optional: plot
    figure
    errorbar(yrs, Trnsprt, Trns_uncert, '-r', 'markersize', 6)
    grid on
    xlabel('Year'); ylabel('Transport (Sv)')
    title('0–500 m Baroclinic Transport into the Nordic Seas')
end