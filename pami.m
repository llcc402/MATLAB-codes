%--------------------------------------------------------------------------
% STEP 0: Init
%--------------------------------------------------------------------------

clear
clc

doc_by_year = [0, 119, 217, 319, 450, 593, 720, 845, 994, 1131, 1262, 1386,...
    1488, 1634, 1786, 1928, 2107, 2289, 2474, 2650, 2833, 3006, 3199, 3380, ...
    3599, 3787, 3979];

% the number of particles
particleN = 1000;

% the number of activated centers
actN = 100;

% the bounds of the symetric-KL
B = 3;

% the concentration parameter of G0
gamma = 5;

% the concentration parameter of G2
alpha = 5;

% the number of Gibbs iterations
maxIter = 500;

%--------------------------------------------------------------------------
% STEP 1: Preprocessing
%--------------------------------------------------------------------------

data = csvread('puredata.csv');

% transform the observations into 12 dimensional vectors
[~, vecs] = sym_cluster(data, 12);

% change the std of the columns of vecs into 1
vecs = vecs ./ repmat(std(vecs), size(vecs, 1), 1);

% one of the element of vecs is extremely large and we modify it into a
% smaller value
ix = vecs < -60;
vecs(ix) = mean(vecs(vecs(:,2) > -60, 2));
vecs(:,2) = vecs(:,2) / std(vecs(:,2));

% init G0 with a DPMM
[z, G0, centers] = dp(vecs, gamma, .1);


% a cell, each value of the cell is corresponding to a year
z_by_year = cell(1, 26);
for i = 1:length(z_by_year)
    z_by_year{i} = z(doc_by_year(i)+1 : doc_by_year(i+1));
end

for iter = 1:maxIter
    
    %--------------------------------------------------------------------------
    % STEP 2: Sample the mixing measure of the first year
    %--------------------------------------------------------------------------

    % the rows are corresponding to the mixing measures of the years
    G2 = zeros(26, length(G0), particleN);

    weight = zeros(1, particleN);

    for i = 1:particleN
        % propose from the prior
        G2(1, :, i) = dpDisrnd(alpha, G0);

        % weighting
        x = G2(1, :, i);
        w = x(z_by_year{1});
        w(w == 0) = 1e-150;
        weight(i) = sum(log(w));
    end

    weight = weight - max(weight);
    weight = exp(weight);

    % normalize the weights
    weight = weight / sum(weight);

    % resample
    [~, ~, ix] = histcounts(rand(1,particleN), [0, cumsum(weight)]);
    G2(1, :, :) = G2(1, :, ix);

    %--------------------------------------------------------------------------
    % STEP 3: Sample the rest of the mixing measures
    %--------------------------------------------------------------------------
    for year = 2:26
        weight = zeros(1, particleN);

        for i = 1:particleN
            % propose from the prior
            G2(year, :, i) = smoothSample(G0, G2(year-1, :, i), B, alpha);

            % weighting
            x = G2(year, :, i);
            w = x(z_by_year{year});
            w(w == 0) = 1e-150;
            weight(i) = sum(log(w));
        end

        weight = weight - max(weight);
        weight = exp(weight);

        % normalize
        weight = weight / sum(weight);

        % resample
        [~, ~, ix] = histcounts(rand(1,particleN), [0, cumsum(weight)]);
        G2(year, :, :) = G2(year, :, ix);
    end

    % use the mean as our sampled mixing measures
    G2_mean = mean(G2, 3);

    %--------------------------------------------------------------------------
    % STEP 4: Sample indicators
    %--------------------------------------------------------------------------

    for year = 1:26

        for j = 1:length(z_by_year{year})

            distance = repmat(vecs(year + j,:), actN, 1) - centers;
            likelihood = - sum(distance .^ 2, 2);

            log_post = log(G2_mean(year, :)) + likelihood';
            log_post = log_post - max(log_post);
            post = exp(log_post);
            post = post / sum(post);

            [~, ~, z_by_year{year}(j)] = histcounts(rand(1), [0, cumsum(post)]);
        end
    end

    %--------------------------------------------------------------------------
    % STEP 5: Sample centers
    %--------------------------------------------------------------------------

    z = cell2mat(z_by_year);
    ix = accumarray(z', 1:length(z), [], @(x){x});
    for i = 1:length(ix)
        if ~isempty(ix{i})
            centers(i,:) = mean(vecs(ix{i}, :));
        end
    end

%     %--------------------------------------------------------------------------
%     % STEP 6: Sample G0
%     %--------------------------------------------------------------------------
% 
%     counts = histcounts(z, 1:actN+1);
%     a = counts + 1;
%     b = [cumsum(counts(2:end), 'reverse'), 0];
%     b = b + gamma;
%     V = betarnd(a, b);
%     G0 = V;
%     V = cumprod(1 - V);
%     G0(2:end) = G0(2:end) .* V(1:end-1);
    
    sample.G0 = G0;
    sample.G2 = G2_mean;
    sample.z = z;
    save([num2str(iter), 'sample.mat'], 'sample')

    fprintf('iteration %d done\n', iter)
end






