%--------------------------------------------------------------------------
% This function is used to sample from a Dirichlet process with discrete
% base measure.
%--------------------------------------------------------------------------
% Input:
%     alpha         a scalar, the concentration parameter.
%     G_0           a discrete measure, the base measure.
% Output:
%     G             a discrete measure having the same atoms with G_0. In
%                   paticular, G is a vector of probabilities here.
%--------------------------------------------------------------------------
% Model:
%     G ~ DP(alpha, G_0)
% Concretely, let G_0 = (p1,...,p_n,...) and G = (c1,...,c_n,...), we have
% the following sampling scheme:
%     c_tilde_n ~ Beta(alpha * p_n, alpha * (1 - sum_{k=1}^n p_n)
%     c_n = c_tilde_n * prod_{k=1}^{n-1} (1 - c_tilde_k)
%
function G = dpDisrnd(alpha, G_0)
if nargin < 2
    G_0 = gem();
end
if nargin < 1
    alpha = 1;
end

param1 = alpha * G_0;
param2 = alpha * (1 - cumsum(G_0));

G_tilde = betarnd(param1, param2);
G = G_tilde;
G_tilde = cumprod(1 - G_tilde);
G(2:end) = G(2:end) .* G_tilde(1:end-1);

% in case 1 - sum(G0) < 0
G(end) = 0;
end