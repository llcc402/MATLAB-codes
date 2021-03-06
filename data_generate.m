% Input:
%     G0       a row vector. The base probability of the DP.
%     B        a scalar. The bound of the kl.
%     K        a scalar. The number of groups.
%     n        a scalar. The number of activated weights in the distros.
%     m        a scalar. The number of observations in each group.
%     alpha    a scalar. The concentration of the DP where Gk sampled from.
%     concent  a scalar. Controls the variance of the noise.
% Output:
%     distro   a matrix of order K * n.
%     data     a matrix of order K * m. 
%--------------------------------------------------------------------------
% Model:
%     data(k,i) ~ distro(k,:), k = 1,...,K, i = 1,...,m 
%     distro    =  (G1; G2; G3; ...; G_K)
%     G_k       ~   G0, k=1, ..., K
%              s.t. KL(G_{k-1}||G_k) <= B
function [distro, distro_noise, data, kl] = data_generate(G0, B, K, n, m, alpha, concent)
if nargin < 7
    concent = 0.03;
end
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

% init
kl = zeros(1, K-1);
distro = zeros(K, n);
distro_noise = distro;

% generate G1
distro(1, :) = dpDisrnd(alpha, G0);
% add noise
distro_noise(1,:) = add_noise(distro(1,:), concent);

for k = 2:K    
    distro(k, :) = smoothSample(G0, distro(k-1, :), B, alpha);
    % add noise
    distro_noise(k,:) = add_noise(distro(k,:), concent);

    kl(k-1) = symKL(distro(k-1,:), distro(k,:));
end

% generate data points
data = zeros(K, m);
for k = 1:K
    data(k, :) = discreternd(m, distro_noise(k, :));
end

end

function noisy_data = add_noise(data, concent)
% data is a row vector
noise = gamrnd(concent, 1, 1, length(data));
noise = noise / sum(noise);
noisy_data = data + noise;
noisy_data = noisy_data / sum(noisy_data);
end