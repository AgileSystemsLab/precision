function [SS] = pointLineSSerr(x, y, beta)
% Quick function to get sum of squared distances between set of points and
% a line defined by slope=beta(1), intercept=beta(2)
SS = abs(beta(1) .* x - y + beta(2)) ./ sqrt(beta(1)^2 + 1);
SS = sum(SS);
end