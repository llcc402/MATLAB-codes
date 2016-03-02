%--------------------------------------------------------------------------
% This function is used to sample G2(1) given the base measure G0 and the 
% precedent G1. 
%--------------------------------------------------------------------------
% Input:
%     G0      a row vector. The base measure
%     G1      a row vector. The precedent measure.
%     B       a scalar. The bound of the symetric kl.
%     alpha   a scalar. The concentration parameter.
% Output:
%     x       a scalar. The weight of the first atom of G1.
% Model:
%     G1 ~ DP(G0)
%     G2 ~ DP(G0)
%     symKL(G1, G2) <= B
%     G2(1) = x

function x = sample_the_first(G0, c, B, alpha)
if nargin < 4
    alpha = 1;
end
if nargin < 3
    B = 1;
end

if c == 0
    error('G1(1) should not be 0!')
end

round = 1e-100; % we lood for roots in the interval [round, 1-round]
pace = 1e-5; % we lood for root in the intervals [i, i + pace]

% func = @(x) x * log(x/c) + (1-x) * log((1-x) / (1-c)) + c * log(c/x) ...
%      + (1-c) * log((1-c) / (1-x)) - B;
% fun_1 = func(round);
% fun_2 = func(1 - round);
% if fun_1 * fun_2 < 0
%     root = fzero(func, [round, 1-round]);
% elseif fun_1 < 0
%     root = [];
% else
%     [x_min, fval] = fminbnd(func, round, 1-round);
%     if fval >= 0
%         error('G2(1) has no solution')
%     else
%         root = zeros(1,2);
%         root(1) = fzero(func, [round, x_min]);
%         root(2) = fzero(func, [x_min, 1-round]);
%     end
% end

z = round:pace:1-round-pace;
y = z .* log(z/c) + (1-z) .* log((1-z) / (1-c)) + c * log(c./z) ...
     + (1-c) * log((1-c) ./ (1-z)) - B < 0;
y = diff(y);
ix = find(y ~= 0);
func = @(x) x * log(x/c) + (1-x) * log((1-x) / (1-c)) + c * log(c/x) ...
     + (1-c) * log((1-c) / (1-x)) - B;
if length(ix) == 2
    root = zeros(1, 2);
    root(1) = fzero(func, [z(ix(1)), z(ix(1)+1)]);
    root(2) = fzero(func, [z(ix(2)), z(ix(2)+1)]);
elseif length(ix) == 1
    root = fzero(func, [z(ix), z(ix+1)]);
else
    root = [];
end

% find the truncating points for the beta distribution
% root = [];
% func = @(x) x * log(x/c) + (1-x) * log((1-x) / (1-c)) + c * log(c/x) ...
%     + (1-c) * log((1-c) / (1-x)) - B;
% for i = round : pace : 1-round-pace
%     try
%         root(end+1) = fzero(func, [i, i + pace]);
%     end
% end

% % plot
% z = round:pace:1-round;
% y = z .* log(z/c) + (1-z) .* log((1-z) / (1-c)) + c * log(c./z) ...
%     + (1-c) * log((1-c) ./ (1-z)) - B;
% plot(z, y, root, zeros(1, length(root)), 'o')

% find the feasible interval
if isempty(root)
    if 1/2 * log(1/2/c) + (1-1/2) * log((1-1/2) / (1-c)) + c * log(c/1*2) ...
       + (1-c) * log((1-c) / (1-1/2)) - B < 0
        interval = [0 1];
    else
        error('something is wrong because i cannot find a root')
    end
elseif length(root) == 1
    if root/2 * log(root/2/c) + (1-root/2) * log((1-root/2) / (1-c)) + ...
       c * log(c/root*2) + (1-c) * log((1-c) / (1-root/2)) - B < 0
        interval = [0 root];
    else
        interval = [root 1];
    end
elseif length(root) == 2
    interval = root;
else
    error('More roots than we expect')
end
    
% sample G2(1)
x = 0;
n = 0;
while x == 0 && n < 99
    x = truncBetarnd(interval(1), interval(2), alpha * G0(1), alpha * (1-G0(1)), 1);
    n = n + 1;
end
if n >= 99
    x = 1e-150;
end

end