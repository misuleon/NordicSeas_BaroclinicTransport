function total = sum_Z(z, prop)
% SUM_Z Integrates a property over depth using the trapezoidal rule.
%
%   total = sum_Z(z, prop)
%
%   Inputs:
%       z    - depth vector (must be column vector, length m)
%       prop - property to integrate (size m x n, where n is number of profiles)
%
%   Output:
%       total - integrated property for each profile (1 x n)
%
%   This routine approximates the vertical integral of 'prop' using the
%   trapezoidal rule between successive levels in 'z'.

    [m, n] = size(prop);
    
    if size(z,1) ~= m
        error('Depth vector z must match number of rows in prop.');
    end

    dz = diff(z);                   % vertical spacing (m-1 x 1)
    DZ = repmat(dz, 1, n);          % replicate for all profiles
    prop_avg = (prop(1:end-1,:) + prop(2:end,:)) / 2;  % trapezoidal average

    total = sum(DZ .* prop_avg, 1);  % integrate along vertical

end