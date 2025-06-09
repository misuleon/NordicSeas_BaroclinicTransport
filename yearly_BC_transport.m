% this code depends upon Nordic_Seas_MOC_analysis.m to prepare the data
% sets used here. The steps are to pull together yearly average dyn. hgt 
% profiles for the NS and RT. These have means and standard deviations.
% The RT-NS difference gives velocity. This is integrated in the vertical
% to yield transport for that year. This is done for all available years.
% We have no choice but to assume the same reference depth for the
% integration for all years. 

% first integrate each dynamic hgt profile in vertical, this is effectively
% getting the potential energy. This way we can get mean and scatter for
% the yearly averages. The yearly difference will give us transport. 
% The dynamic height profiles are dyn_NS and dyn_RT but first reference to
% 500 m, i.e. lvl 25:
dyn_NS_ = dyn_NS(25,:) - dyn_NS;
dyn_RT_ = dyn_RT(25,:) - dyn_RT;
PE_NS = sum_Z(dpth,dyn_NS_);
PE_RT = sum_Z(dpth,dyn_RT_);
HT_NS = sum_Z(dpth,temp_NS);
HT_RT = sum_Z(dpth,temp_RT);

% eliminate two RT stations in 1949. They are excessively warm as if
% profiles in a meddy
chk = 0;
if chk == 0
    lnrt = length(yrRT);
    yrRT = yrRT([1:16 19:lnrt]);
    PE_RT = PE_RT([1:16 19:lnrt]);
    chk = 1;
end
%%%%%% define period of interest %%%%%%
yrSt = 1900; yrEn = 2023;

N_yrs = yrEn - yrSt + 1;
for i = 1:N_yrs
    iy = yrSt + i - 1;
% first NS
    k = find(yrNS == iy);
    mean_PE_NS(i) = NaN; std_PE_NS(i) = NaN; N_NS(i) = NaN;
    mean_PE_NS(i) = mean(PE_NS(k));    % PE at sfc rel 475 m
    std_PE_NS(i) = std(PE_NS(k));
    N_NS(i) = length(k);
    date_year_NS(i) = datenum(iy,7,1);
    
% RT
    j = find(yrRT == iy);
    mean_PE_RT(i) = NaN; std_PE_RT(i) = NaN; N_RT(i) = NaN;
    mean_PE_RT(i) = mean(PE_RT(j));    % PE at sfc rel 475 m
    std_PE_RT(i) = std(PE_RT(j));
    N_RT(i) = length(j);
    date_year_RT(i) = datenum(iy,7,1);
end

% move NS 1936 to 1938 to coincide with RT
    mean_PE_NS(39) = mean_PE_NS(37);
    std_PE_NS(39) = std(PE_NS(37));
    
for i = 1:N_yrs
    iy = yrSt + i - 1;
% now difference and uncertainty
    Trnsprt(i) = NaN; Trns_uncert(i) = NaN;
    Trnsprt(i) = mean_PE_RT(i) - mean_PE_NS(i);
    E = sqrt((std_PE_NS(i).^2 + std_PE_RT(i).^2));
    Trns_uncert(i) = E/sqrt(N_NS(i) + N_RT(i));
end

Trnsprt = Trnsprt*10/eff/1e6;
Trns_uncert = Trns_uncert*10/eff/1e6;
yrs = yrSt:yrEn;

figure(35)
subplot(3,1,1)
errorbar(yrs,Trnsprt,Trns_uncert,'-r','markersize',10)
% plot(yrs,Trnsprt,'-xk','markersize',5)
axis([yrSt yrEn 4 8])
grid on
% xlabel('year','fontsize',16)
ylabel('Sv','fontsize',16)
% title('0-500 m baroclinic transport into Nordic Seas','fontsize',16)
hold on


if 0
    J = find(~isnan(Trnsprt) & yrs > 1949);
    P = polyfit(yrs(J),Trnsprt(J),1)
    yy = [1950 2022];
    tt = polyval(P,yy);
    plot(yy,tt,'k','linewidth',2)
    yy = [1902 1950];
    tt = polyval(P,yy);
    plot(yy,tt,'--k','linewidth',1)
    axis([1900 2022 4 8])
end

if 0
    J = find(~isnan(Trnsprt) & yrs >= 1949);
    yrsJ = yrs(J);
    mid_tid = (yrsJ(1) + yrsJ(end))/2;
    [A,Aerr,xp,yp,yperr ] = linregress(yrs(J), Trnsprt(J));
    % xp = xp + mid_tid;
    figure(35)
    subplot(3,1,1)
    disp('get here?')
    h = plot(xp,yp,'k',xp,yp+yperr,'k--',xp,yp-yperr,'k--','linewidth',1);
    x_x = [1900 1950];
    y_x = [A(1) + 1900*A(2) A(1) + 1950*A(2)]
    plot(x_x,y_x,'k--','linewidth',1)
end

% now get heat variance in RT
for i = 1:N_yrs
    iy = yrSt + i - 1;
    j = find(yrRT == iy);
    mean_HT_RT(i) = NaN; std_HT_RT(i) = NaN; N_RT(i) = NaN;
    mean_HT_RT(i) = mean(HT_RT(j));    % mean temp  0-500 m
    std_HT_RT(i) = std(HT_RT(j));
    N_RT(i) = length(j);
    date_year_RT(i) = datenum(iy,7,1);
    std_err_RT(i) = std(HT_RT(j))/length(j);
end

% now plot normalized variability 51 = 1950

avg_HT_RT = nanmean(mean_HT_RT(51:end));
std_HT_RT = nanstd(mean_HT_RT(51:end));
disp('RT mean and std:')
[avg_HT_RT/500 std_HT_RT/500]

norm_RT = (mean_HT_RT - avg_HT_RT)/std_HT_RT;
std_err_RT = std_err_RT/std_HT_RT;
subplot(3,1,2)
errorbar(yrs,norm_RT,std_err_RT,'-r','markersize',10)
% plot(yrs,norm_RT,'-xr','markersize',5)
grid on
%axis([yrSt yrEn -3 3])
axis([1900 2023 -3 3])
h = text(1910,1,'9.91 ± 0.32°C');
set(h,'fontsize',14,'color','r')

% xlabel('year','fontsize',16)
ylabel('temperature','fontsize',16)
% title(' normalized temperature variability in Rockall Trough','fontsize',16)

% now get heat variance in NS
for i = 1:N_yrs
    iy = yrSt + i - 1;
    j = find(yrNS == iy);
    mean_HT_NS(i) = NaN; std_HT_NS(i) = NaN; N_NS(i) = NaN;
    mean_HT_NS(i) = mean(HT_NS(j));    % mean temp  0-500 m
    std_HT_NS(i) = std(HT_NS(j));
    N_NS(i) = length(j);
    date_year_NS(i) = datenum(iy,7,1);
    std_err_NS(i) = std(HT_NS(j))/length(j);

end

avg_HT_NS = nanmean(mean_HT_NS(51:end));
std_HT_NS = nanstd(mean_HT_NS(51:end));
disp('NS mean and std:')
[avg_HT_NS/500 std_HT_NS/500]

norm_NS = (mean_HT_NS - avg_HT_NS)/std_HT_NS;
std_err_NS = std_err_NS/std_HT_NS;

% move NS 1936 to 1938 to coincide with RT
norm_NS(39) = norm_NS(37);
norm_NS(37) = NaN;

subplot(3,1,3)
errorbar(yrs,norm_NS,std_err_NS,'-b','markersize',10)
% this line to inject this plot into plot_transport_and_trend.m
if 1
    errorbar(yrs+0.5,norm_NS,std_err_NS,'-b','markersize',10)
end
% plot(yrs,norm_NS,'-xb','markersize',5)
grid on
%axis([yrSt yrEn -3 3])
axis([1900 2023 -3 3])
h = text(1910,1,'2.65 ± 0.35°C');
set(h,'fontsize',14,'color','b')
% xlabel('year','fontsize',16)
ylabel('temperature','fontsize',16)
% title(' normalized temperature variability in Norwegian Sea','fontsize',16)

% readme = 'RT = Tockall Trough NS = Norwegian Sea HT = heat content (0-500 m vert avg temp)';
% save save_to_leon_rev_2 readme yrs Trnsprt Trns_uncert mean_HT_RT mean_HT_NS


% get correlation between Trnsprt and heat content in RT and NS
[yrs' norm_RT' norm_NS' Trnsprt'];

% limit to years 1949 to 2022

J = 50:122;
norm_RT_ = norm_RT(J); 
norm_NS_ = norm_NS(J); 
Trnsprt_ = Trnsprt(J);
K = ~isnan(norm_RT_);
norm_RT_ = norm_RT_(K); 
norm_RT_ = norm_RT_/std(norm_RT_);
norm_NS_ = norm_NS_(K); 
Trnsprt_ = Trnsprt_(K);
norm_Trn = (Trnsprt_ - mean(Trnsprt_))/std(Trnsprt_);

Cov_Trn_RT = cov(norm_Trn,norm_RT_)

Cov_Trn_NS = cov(norm_Trn,norm_NS_)

Cov_RT_NS  = cov(norm_RT_/std(norm_RT_),norm_NS_/std(norm_NS_))

