% Input:
%     data        a matrix of order m by n, where m is the number of groups 
%                 and n is the number of observations in each group.
%     gibsN       a scalar. The number of gibbs sampling iteration.
%     particleN   a scalar. The number of particles.
%     actN        a scalar. The number of activated weights in G0.
%     alpha       a scalar. The concentration parameter of G1,...,G_N.
%     B           a scalar. The bound of the symKL.   
%     gamma       a scalar. The concentration parameter of G0
% Output:
%     distro      a matrix of order m * actN, where m is number of groups. 
%     G0          a row vector of length actN.

function [distro, G0] = synthc_particle(data, gibsN, particleN, actN, alpha, B, gamma)
if nargin < 7
    gamma = 5;
end
if nargin < 6
    B = 1;
end
if nargin < 5
    alpha = 1;
end
if nargin < 4
    actN = 100;
end
if nargin < 3
    particleN = 1000;
end
if nargin < 2
    gibsN = 10;
end

%--------------------------------------------------------------------------
% STEP 1: Init
%--------------------------------------------------------------------------
counts = histcounts(data);
if length(counts) < actN
    counts = [counts, zeros(1, actN - length(counts))];
end
G0 = dirichletrnd(gamma + counts);

%--------------------------------------------------------------------------
% STEP 2: Gibbs sampling
%--------------------------------------------------------------------------
for gibsIter = 1:gibsN
    
    % particle filtering method sampling distro
    distro = particle(data, G0, alpha, particleN, B);
    
    % sample G0
    G0 = dirichletrnd(gamma + counts);
    G0 = G0 / sum(G0);
end

end

function distro = particle(data, G0, alpha, particleN, B)

dist = zeros(size(data, 1), length(G0), particleN);

%--------------------------------------------------------------------------
% STEP 1: sample G1
%--------------------------------------------------------------------------

% sampling
% note in an array x(:,1,1) is a column, x(1,:,1) is a row, x(1,1,:) is
% neither a row nor a column
for i = 1:particleN
    dist(1, :, i) = dpDisrnd(alpha, G0);
end

% weighting
log_weight = zeros(1, particleN);
for i = 1:particleN
    log_weight(i) = get_weight(dist(1, :, i), data(1, :));
end
log_weight = log_weight - mean(log_weight(log_weight ~= -Inf));
weight = exp(log_weight);

% normalize
if sum(weight) == 0
    error('the weights of the first distro are all 0')
end
weight = weight / sum(weight);

% resampling
[~, ~, ix] = histcounts(rand(1, particleN), [0, cumsum(weight)]);
dist(1,:,:) = dist(1, :, ix);

%--------------------------------------------------------------------------
% STEP 2: sample G_{2:end}
%--------------------------------------------------------------------------
for jj = 2:size(data,1)
    % sampling
    for i = 1:particleN
        dist(jj, :, i) = smoothSample(G0, dist(jj-1, :, i), B, alpha);
    end
    
    % weighting
    for i = 1:particleN
        log_weight(i) = get_weight(dist(jj, :, i), data(jj, :));
    end
    if sum(exp(log_weight)) == 0
        g = dist(jj-1, :, i);
        log_weight = log_weight - median(log_weight(log_weight ~= -Inf));
    end
    weight = exp(log_weight);
    
    if sum(weight) == Inf
        error('one of the weight is +Inf')
    end
    if sum(weight) == 0
        error(['the weights of distro ',num2str(jj), ' are all 0'])
    end
    if sum(isnan(weight))
        error('one of the weight is nan')
    end
    weight = weight / sum(weight);
    
    % resample
    [~, ~, ix] = histcounts(rand(1, particleN), [0 cumsum(weight)]);
    dist(jj, :, :) = dist(jj, :, ix);    
end

%--------------------------------------------------------------------------
% STEP 3: mean as sample of distro
%--------------------------------------------------------------------------
distro = mean(dist, 3);

end

function log_weight = get_weight(d, data)
w = d(data);
w(w == 0) = 1e-150;
log_weight = sum(log(w));

end

function x = dirichletrnd(gamma)
% the input alpha is a row vector
% the output x is a row vector of the same size of gamma
y = gamrnd(gamma, 1);
x = y ./ sum(y);

end




















