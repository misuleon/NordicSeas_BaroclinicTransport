function [A, amp] = deseason(tim,var,fmt)
% routine deseason.m first gets monthly means and then removes these.
% fmt = 1: tim is matlab time (days)
% fmt = 2: tim is in years.decimals
% first get month of variable
% tim =  datenum(1950,1,1) + [1:1:365] - 1;
% var = 0.5*sin((tim-tim(1))/365*2*pi) + 6;
% fmt = 1;
fmt;
if fmt == 1
    [yr,mo,dy] = datevec(tim);
elseif fmt == 2
    mo = fix(rem(tim,1)*12) + 1;
end

%first remove mean (add back later)
N = ~isnan(var);
mean_var = mean(var(N));
var_rest = var - mean_var;
% now get means for each month:
for i = 1:12
    J = mo == i;
    sea_amp(i) = nanmean(var_rest(J));
end

amp = std(sea_amp);
% add back mean and remove monthly signal
A = N.*(var - sea_amp(mo));

% plot(tim,var,'k',tim,A,'*r')
