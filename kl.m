function r = kl(p,q)
if nargin < 1
    p = rand(1,100);
    p = p / sum(p);
    q = rand(1,100);
    q = q / sum(q);
end
if nargin < 2
    if min(size(p)) > 1
        error('The arguments should be vectors!');
    else
        q = rand(size(p));
        q = q / sum(q);
    end
end
    
if min(size(p)) > 1 || min(size(q)) > 1
    error('The arguments should be vectors!');
end

if size(p, 1) > 1
    p = p';
end
if size(q, 1) > 1
    q = q';
end

if sum(q) > 1
    q = q / sum(q);
end
if sum(p) > 1
    p = p / sum(p);
end

% epsilon is used to discard small numbers, in case the kl goes to inf
epsilon = 1e-10;

p(p < epsilon) = epsilon;
q(q < epsilon) = epsilon;
r = log(p ./ q);
r = p * r';

% if (sum(p < epsilon) == 0) && (sum(q < epsilon) == 0)
%     r = log(p./q);
%     r = p * r';
% else
%     ix_p = find(p > epsilon);
%     ix_q = find(q > epsilon);
%     ix = intersect(ix_p, ix_q);
%     
%     % if the sum of the intersection of atoms taking positive probabilities 
%     % of p and q less than 0.95, we believe that the kl of p and q should
%     % be infinite.
%     if sum(p(ix)) < 0.9 || sum(q(ix)) < 0.9
%         r = inf;
%     elseif sum(p(ix)) < 0.95 || sum(q(ix)) < 0.95
%         warning('The result is not precise because some atoms are neglected');
%         r = log(p(ix) ./ q(ix));
%         r = p(ix) * r';
%     else
%         r = log(p(ix) ./ q(ix));
%         r = p(ix) * r';
%     end
% end

end