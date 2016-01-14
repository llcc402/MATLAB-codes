function distro_trad_mean = traditional_particle(data, particleN, actN, gamma, alpha, gibsN)
if nargin < 6
    gibsN = 1;
end
if nargin < 5
    alpha = 1;
end
if nargin < 4
    gamma = 0;
end

% init
distro_trad = zeros(size(data,1), actN);
distro_trad_mean = zeros(size(distro_trad));

for iter = 1:gibsN
    % generate G0
    counts = histcounts(data, 1:actN+1);
    G0 = dirichletrnd(counts + gamma);

    % generate G(1:end)
    for i = 1:size(data,1)
        distro_trad(i,:) = particle(G0, data(i,:), particleN, alpha, actN);
    end
    distro_trad_mean = distro_trad_mean + distro_trad;
end
distro_trad_mean = distro_trad_mean / gibsN;

end

%%
function dist = particle(G0, data_vec, particleN, alpha, actN)
% init
dist = zeros(particleN, actN);

% sampling
for i = 1:size(dist, 1)
    dist(i, :) = dpDisrnd(alpha, G0);
end

% weighting
log_weight = zeros(1, particleN);
for i = 1:particleN
    log_weight(i) = get_weight(dist(i,:), data_vec);
end
if sum(exp(log_weight)) == 0
    log_weight = log_weight - max(log_weight(log_weight ~= -Inf));
end
weight = exp(log_weight);
weight = weight / sum(weight);

% resampling
[~, ~, ix] = histcounts(rand(1, particleN), [0, cumsum(weight)]);
dist = dist(ix,:);

% mean as posst
dist = mean(dist);

end

%%
function r = dirichletrnd(gamma)
r = gamrnd(gamma, 1);
r = r / sum(r);
end

%%
function w = get_weight(prob, data)
w = prob(data);
w(w == 0) = 1e-200;
w = log(w);
w = sum(w);
end
