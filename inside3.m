function Bnan = inside3(xx, yy, xb, yb)
% INSIDE3  Returns a mask for points inside a polygon defined by (xb, yb).
%
%   Bnan = inside3(xx, yy, xb, yb)
%
%   Inputs:
%       xx, yy - arrays of coordinates (can be vectors or matrices)
%       xb, yb - vectors defining the polygon boundary
%
%   Output:
%       Bnan   - mask of same size as xx and yy:
%                  1 for points inside polygon
%                  NaN for points outside
%
%   Notes:
%     - Uses angle summation method to determine point-in-polygon.
%

    nb = length(xb);
    x = xx(:);
    y = yy(:);
    np = numel(x);
    Bnan = nan(size(xx));

    for p = 1:np
        if x(p) >= min(xb) && x(p) <= max(xb) && ...
           y(p) >= min(yb) && y(p) <= max(yb)

            alpha = zeros(nb-1,1);

            for b = 2:nb
                % Triangle sides for Law of Cosines
                a2 = (xb(b)   - xb(b-1))^2 + (yb(b)   - yb(b-1))^2;
                b2 = (xb(b-1) - x(p))^2   + (yb(b-1) - y(p))^2;
                c2 = (xb(b)   - x(p))^2   + (yb(b)   - y(p))^2;

                if b2 ~= 0 && c2 ~= 0
                    alpha(b-1) = acos(0.5 * (b2 + c2 - a2) / sqrt(b2 * c2));

                    % Adjust angle sign based on orientation
                    phi = atan2(yb(b) - yb(b-1), xb(b) - xb(b-1));
                    yp = -(xb(b) - x(p)) * sin(phi) + (yb(b) - y(p)) * cos(phi);
                    if yp < 0
                        alpha(b-1) = -alpha(b-1);
                    end
                end
            end

            % If sum of angles ≈ 2π (or -2π), point is inside
            if abs(sum(alpha)) > 3
                Bnan(p) = 1;
            end
        end
    end
end