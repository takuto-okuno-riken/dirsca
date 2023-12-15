% get nuisance time-series of Linear, quadratic drift and constant term
% returns 2 linear, 4 quadratic drift, 1 constant time-series.
% input:
%  n            time-series length

function Xn = getNuisancePolynomial(n)
    Xn = nan(n,7);
    Xn(:,1) = 0:1/(n-1):1;
    Xn(:,2) = 1:-1/(n-1):0;
    Xn(:,3) = Xn(:,1) .* Xn(:,1);
    Xn(:,4) = Xn(end:-1:1,3);
    Xn(:,5) = 1 - Xn(:,3);
    Xn(:,6) = 1 - Xn(:,4);
    Xn(:,7) = 1;
end
