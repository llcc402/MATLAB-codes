%--------------------------------------------------------------------------
% The symetric KL of discrete probability distributions p and q.
% Input:
%     p     a row vector.
%     q     a row vector. Note p and q should have the same length
% Output:
%     r     a scalar. The symetric KL of p and q.
% Model:
%     r = KL(p,q) + KL(q,p)

function r = symKL(p, q)
if length(p) ~= length(q)
    error('The length of p and q should equal')
end

if size(p, 1) > 1
    error('p should be a row vector')
end
if size(q, 1) > 1
    error('q should be a row vector')
end

% ix_p = find(p > 1e-50);
% ix_q = find(q > 1e-50);
% ix = intersect(ix_p, ix_q);
% 
% if sum(p(ix)) < 0.9 || sum(q(ix)) < 0.9
%     warning('Some of the probabilities are neglected')
% end
% 
% p = p(ix);
% q = q(ix);
% 
% r = p * (log(p./q))' + q * (log(q./p))';
epsilon = 1e-5;
p(p < epsilon) = epsilon;
q(q < epsilon) = epsilon;

r = p * (log(p./q))' + q * (log(q./p))';

end