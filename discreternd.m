% Note the output x is a row vector
function x = discreternd(n, probs)
if size(probs, 1) > 1 % change probs from column vector to row vector
    probs = probs';
end
r = rand(1, n);
probs = [0 cumsum(probs)];

[~, ~, x] = histcounts(r, probs);
end