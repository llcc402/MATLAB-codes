function w = gem(n, alpha)
if nargin < 2
    alpha = 1;
end
if nargin < 1
    n = 100;
end

v = betarnd(1, alpha * ones(1,n));
w = v;
v = cumprod(1 - v);
w(2:end) = w(2:end) .* v(1:end-1);
end