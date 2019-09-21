function [ diff ] = projdiff( tau,alpha,projpoints,littlex )
xydiff = gsubtract(littlex,projpoints');
diff = sum( alpha .* ((tau^2)/2 .* sign(xydiff) .* exp(-tau .* abs(xydiff) ./ sqrt(2)) .* sin(tau .* abs(xydiff) ./ sqrt(2))));
end

