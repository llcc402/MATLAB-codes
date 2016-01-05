%--------------------------------------------------------------------------
% We try to sample G_2 from a DP with base measure G_0 that 
%                    kl(G_1||G_2) < B
% Note: sometimes, the sampled G2 is not a probability distribution in the
% sense that the sum of G2 is smaller than 1. In this case, we should
% discard this sample.
%--------------------------------------------------------------------------
% Input:
%     G0        a row vector. The base measure of G1 and G2.
%     G1        a row vector. 
%     B         a scalar. The bound.
%     alpha     a scalar. The concentration of G1 and G2.
%     pace      a scalar. We test roots in the intervals [i, i + pace]
% Output:
%     G2        a row vector.
%--------------------------------------------------------------------------
function G2 = smoothSample(G0, G1, B, alpha)
if nargin < 4
    alpha = 1;
end
if nargin < 3
    B = 1;
end
if nargin < 2
    G1 = gem();
end
if nargin < 1
    G0 = gem();
end

G1 = G1 / sum(G1);
ix = find(G1 > 1e-100);
G2 = zeros(1, length(G1));

if length(ix) == 1
    error('length(G1 > 1e-100) = 1')
end

G2(ix(1)) = sample_the_first(G0, G1(ix(1)), B, alpha);
for k = 2:length(ix)
    G2(ix(k)) = sample_the_rest(G0, G1, G2, ix(k), B, alpha);
end

end




















































