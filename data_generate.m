% Input:
%     G0       a row vector. The base probability of the DP.
%     B        a scalar. The bound of the kl.
%     K        a scalar. The number of groups.
%     n        a scalar. The number of activated weights in the distros.
%     m        a scalar. The number of observations in each group.
%     alpha    a scalar. The concentration of the DP where Gk sampled from.
% Output:
%     distro   a matrix of order K * n.
%     data     a matrix of order K * m. 
%--------------------------------------------------------------------------
% Model:
%     data(k,i) ~ distro(k,:), k = 1,...,K, i = 1,...,m 
%     distro    =  (G1; G2; G3; ...; G_K)
%     G_k       ~   G0, k=1, ..., K
%              s.t. KL(G_{k-1}||G_k) <= B
function [distro, data] = data_generate(G0, B, K, n, m, alpha)
if nargin < 6
    alpha = 1;
end
if nargin < 5
    m = 50;
end
if nargin < 4
    n = 100;
end
if nargin < 3
    K = 20;
end
if nargin < 2
    B = 1;
end
if nargin < 1
    G0 = gem(n, 5);
end

% generate distributions with base G0
distro = zeros(K, n);
distro(1, :) = dpDisrnd(alpha, G0);
for k = 2:K
    s = 0;
    n = 0;
    while s < 0.99 && n < 50
        distro(k, :) = smoothSample(G0, distro(k-1, :), B, alpha);
        s = sum(distro(k, :));
        n = n + 1;
    end
end

% generate data points
data = zeros(K, m);
for k = 1:K
    data(k, :) = discreternd(m, distro(k, :));
end

end

% G0 = gem(100, 5);
% data = zeros(5, 100);
% data(1, :) = dpDisrnd(1, G0);
% for k = 2:5
%     data(k,:) = smoothSample(G0, G1);
%     kl(data(k-1,:), data(k, :))
% end
