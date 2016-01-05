%--------------------------------------------------------------------------
% This function is used to sample one of G2(2:end). 
%--------------------------------------------------------------------------
% Input:
%     G0      a row vector. The base measure of G1 and G2
%     G1      a row vector. The precedent
%     G2      a row vector. The first k-1 probabilities are known 
%     k       a scalar. The k-th atom of G2
%     B       a scalar. The bound of the symetric kl.
%     alpha   a scalar. The concentration parameter.
% Output:
%     x       a scalar. A probability of some atom in {2, 3, ...}
% Model:
%     G1 ~ DP(alpha, G0)
%     G2 ~ DP(alpha, G0)
%     symKL(G1, G2) <= B
%     G2(k) = x, k \neq 1

function x = sample_the_rest(G0, G1, G2, k, B, alpha)
if nargin < 6
    alpha = 1;
end
if nargin < 5
    B = 1;
end

A1 = sum(G1(1:k-1)) * log(sum(G1(1:k-1))/sum(G2(1:k-1))) ...
     + sum(G2(1:k-1)) * log(sum(G2(1:k-1))/sum(G1(1:k-1)));
c = G1(k);
A2 = 1 - sum(G1(1:k));
A3 = 1 - sum(G2(1:k-1));

if sum(G2(1:k-1)) == 0
    error('sum(G2(1:k-1)) == 0')
end
if sum(G1(1:k-1)) == 0
    error('sum(G1(1:k-1)) == 0')
end

if A2 <= 0 % since all G1(k+1:end) == 0 and G2(k+1:end) should be 0 concequently
    x = 1 - sum(G2(1:k-1));
else
    round = 1e-10; % we look for roots in the interval [round, A3-round] 
    pace = 1e-3; % we look for roots  in the intervals [i, i+pace]

%   % find truncating points
%     func = @(x) A1 + x * log(x/c) + c * log(c/x) + A2 * log(A2 / (A3 - x)) ...
%           + (A3 - x) * log((A3 - x) / A2) - B;
%     fun_1 = func(round);
%     fun_2 = func(A3-round);
%     if fun_1 * fun_2 < 0
%         root = fzero(func, [round, A3-round]);
%         fprintf(['there is only one root ', num2str(root), '\n'])
%         fprintf(['A3 = ', num2str(A3), '\n'])
%     elseif fun_1 < 0
%         root = [];
%     else
%         [x_min, fval] = fminbnd(func, round, A3 - round);
%         if fval > 0
%             root = [];
%         else
%             root = zeros(1, 2);
%             root(1) = fzero(func, [round, x_min]);
%             root(2) = fzero(func, [x_min, A3-round]);
%             fprintf(['the two roots are', num2str(root), '\n'])
%         end
%     end
%     
    z = round:pace:A3-round;
    y = A1 + z .* log(z/c) + c .* log(c./z) + A2 .* log(A2 ./ (A3 - z)) ...
       + (A3 - z) .* log((A3 - z) / A2) - B;
    yy = y < 0;
    yy = diff(yy);
    ix = find(yy ~= 0);
    func = @(x) A1 + x * log(x/c) + c * log(c/x) + A2 * log(A2 / (A3 - x)) ...
          + (A3 - x) * log((A3 - x) / A2) - B;
    if length(ix) == 2
        root = zeros(1, 2);
        root(1) = fzero(func, [z(ix(1)), z(ix(1)+1)]);
        root(2) = fzero(func, [z(ix(2)), z(ix(2)+1)]);
    elseif length(ix) == 1
        root = fzero(func, [z(ix), z(ix+1)]);
    else
        root = [];
    end


%     % plot
%     z = round:pace:A3;
%     y = A1 + z .* log(z/c) + c .* log(c./z) + A2 .* log(A2 ./ (A3 - z)) ...
%         + (A3 - z) .* log((A3 - z) / A2) - B;
%     plot(z, y, root, zeros(1, length(root)), 'o')

    x_not_sampled = true;
    % find feasible intervals
    if isempty(root)
        if A1 + A3/2 * log(A3/2/c) + c * log(c/A3*2) + A2 * log(A2 / (A3/2)) ...
           + (A3/2) * log((A3/2) / A2) - B < 0
            interval = [0 A3];
        else
            x = 0;
%             if G1(k) > 1e-5
%                 fprintf([num2str(y(ix)), ' ', num2str(y(ix+1)), '\n'])
%                 fprintf(['A1 = ', num2str(A1), '\n'])
%                 fprintf(['A2 = ', num2str(A2), '\n'])
%                 fprintf(['A3 = ', num2str(A3), '\n'])
%                 fprintf(['c = ', num2str(c), '\n'])
%                 fprintf(['z1 = ', num2str(z(ix)), '\n'])
%                 fprintf(['z2 = ', num2str(z(ix+1)), '\n'])
%             end
            x_not_sampled = false;
        end
    elseif length(root) == 1
        if A1 + root/2 * log(root/2/c) + c * log(c/root*2) + A2 * log(A2 / (A3 - root/2)) ...
           + (A3 - root/2) * log((A3 - root/2) / A2) - B < 0
            interval = [0 root];
        else
            interval = [root A3];
        end
    elseif length(root) == 2
        interval = root;
    else
        error('More roots than we expected')
    end

    if x_not_sampled
        % change the interval to G2(k)_tilde
        interval = interval / A3;
        % sample G2(k)_tilde
        x = truncBetarnd(interval(1), interval(2), alpha * G0(k), alpha * (1 - sum(G0(1:k))), 1);
        % compute G2(k)
        x = x * A3;
    end
end

end
