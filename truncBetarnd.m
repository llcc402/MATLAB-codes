%--------------------------------------------------------------------------
% Sampling from a truncated Beta distribution.
%--------------------------------------------------------------------------
% Input:
%     a        a scalar, the smaller truncating point
%     b        a scalar, the larger truncating point
%     alpha    a scalar, the first parameter of the beta distribution
%     beta     a scalar, the second parameter of the beta distribbution
%     n        a scalar, the size of the sample
% Output:
%     x        a scalar, sampled from the truncated beta distribution
%--------------------------------------------------------------------------
% Model:
%     x ～ Beta(alpha, beta) * I(a < x < b),
% where 0 < a < b < 1.
%--------------------------------------------------------------------------
% Method:
% Use slice sampler with the joint 
%     p(x,u|alpha,beta) propto x^(alpha-1) I(0 < u < (1-x)^(beta-1), a < x < b)
function x = truncBetarnd(a,b,alpha,beta, n)
if nargin < 5
    n = 5000;
end
if nargin < 4
    beta = 1/3;
end
if nargin < 3
    alpha = 2;
end
if nargin < 2
    b = 0.8;
end
if nargin < 1
    a = 0.2;
end
if beta == 1 || alpha == 1
    error('Please use inverse transform directly!');
end
if ~(0 <= a && a < b && b <= 1)
    fprintf(['a = ', num2str(a), '\n'])
    fprintf(['b = ', num2str(b), '\n'])
    error('Your truncating points are illegal!');
end

if n > 1
    x = zeros(1, n);
    x(1) = rand(1) * (b - a) + a;
    if beta > 1
        for i = 2:n
            u = rand(1) * (1 - x(i-1))^(beta - 1);
            lower= a;
            upper = min(1 - u^(1 / (beta - 1)), b);
            t = rand(1);
            x(i) = nthroot(t * (upper^alpha - lower^alpha) + lower^alpha, alpha);
        end
    elseif beta < 1 && 0 < beta
        for i = 2:n
            u = rand(1) * (1 - x(i-1))^(beta - 1);
            upper = b;
            lower = max(1 - u^(1 / (beta - 1)), a);
            t = rand(1);
            x(i) = nthroot(t * (upper^alpha - lower^alpha) + lower^alpha, alpha);
        end
    end
elseif n == 1
    m = 10;
    x = zeros(1, m);
    x(1) = rand(1) * (b - a) + a;
    if beta > 1
        for i = 2:m
            u = rand(1) * (1 - x(i-1))^(beta - 1);
            lower= a;
            upper = min(1 - u^(1 / (beta - 1)), b);
            t = rand(1);
            x(i) = nthroot(t * (upper^alpha - lower^alpha) + lower^alpha, alpha);
        end
    elseif beta < 1 && 0 < beta
        for i = 2:m
            u = rand(1) * (1 - x(i-1))^(beta - 1);
            upper = b;
            lower = max(1 - u^(1 / (beta - 1)), a);
            t = rand(1);
            x(i) = nthroot(t * (upper^alpha - lower^alpha) + lower^alpha, alpha);
        end
    end
    x = x(m);
end
    
% hist(x, 50)
% title('The histogram of the desired truncated Beta distribution.')
% 
% figure(2)
% s = betarnd(alpha, ones(1,n) * beta);
% s = s(s < b);
% s = s(s > a);
% hist(s, 50)
% title('The histogram of the Beta distribution with the same parameters.')
% 
% p = a:(b-a)/100:b;
% theoretic = (betacdf(p, alpha, beta) - betacdf(a, alpha, beta)) ...
%     / (betacdf(b, alpha, beta) - betacdf(a, alpha, beta));
% exp = zeros(size(p));
% for i = 1:length(p)
%     exp(i) = sum(x < p(i)) / length(x);
% end
% figure(3)
% plot(theoretic, theoretic, theoretic, exp, 'o')
% xlabel('Theoreical cumulative probability')
% ylabel('Experimental cumulative probabilities')

end