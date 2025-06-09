function [Bnan] = inside3(xx,yy,xb,yb)
% function [Bnan] = inside3(xx,yy,xb,yb)
% function to blank a region outside a defined polygon whose vertices are
% defined by the (xb,yb) coordinates;
%
% This is a completely new method from Julio's earlier
% scheme          CNF 3/11/98

nb = length(xb);
x = xx(:); y = yy(:); % Vectorize the positions
np = length(x);
Bnan = nan*ones(size(xx));

for p=1:np
  if (x(p) <= max(xb) & x(p) >= min(xb) & y(p) <= max(yb) & y(p) >= min(yb)) 
    for b = 2:nb
	% Use Law of Cosines to determine angle between
	% the point and two adjacent boundary points
      a2 = (xb(b)-xb(b-1))^2 + (yb(b)-yb(b-1))^2;
      b2 = (xb(b-1)-x(p))^2 + (yb(b-1)-y(p))^2;
      c2 = (xb(b)-x(p))^2 + (yb(b)-y(p))^2;
      if (b2 ~= 0 & c2 ~= 0)
        alpha(b-1) = acos(0.5*(b2 + c2 - a2)/sqrt(b2*c2));
	% Determine the sign of the angle
        phi = atan2(yb(b)-yb(b-1), xb(b)-xb(b-1));
        yp = -(xb(b)-x(p))*sin(phi) + (yb(b)-y(p))*cos(phi);
        if (yp < 0), alpha(b-1) = -alpha(b-1);, end;
     else
        alpha(b-1) = 0;
      end 
    end
	% If the sum of the angles == 2*pi then the point is
	% within the boundary point
    %fprintf('%f\n',(abs(sum(alpha))))
    if (abs(sum(alpha)) > 3)
      Bnan(p) = 1;
    end
  end
end

