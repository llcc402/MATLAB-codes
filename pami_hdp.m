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

% the mixing measures
G2 = zeros(26, length(G0));

% a cell, each value of the cell is corresponding to a year
z_by_year = cell(1, 26);
for i = 1:length(z_by_year)
    z_by_year{i} = z(doc_by_year(i)+1 : doc_by_year(i+1));
end



for iter = 1:maxIter
    %----------------------------------------------------------------------
    % STEP 2: Sampling mixing measure for each year
    %----------------------------------------------------------------------
    for j = 1:26
        counts = histcounts(z_by_year{j}, 1:1+actN);
        a = counts + alpha * G0;
        b = [cumsum(a(2:end), 'reverse'), 0];
        V = betarnd(a, b);
        G2(j,:) = V;
        V = cumprod(1 - V);
        G2(j,2:end) = G2(j,2:end) .* V(1:end-1);
    end
    
    %----------------------------------------------------------------------
    % STEP 3: Sampling z
    %----------------------------------------------------------------------
    for year = 1:26
        
        for j = 1:length(z_by_year{year})
            distance = repmat(vecs(year + j,:), actN, 1) - centers;
            likelihood = - sum(distance .^ 2, 2);
            
            log_post = log(G2(year, :)) + likelihood';
            log_post = log_post - max(log_post);
            post = exp(log_post);
            post = post / sum(post);
            
            [~, ~, z_by_year{year}(j)] = histcounts(rand(1), [0, cumsum(post)]);
        end
    end
    
    %--------------------------------------------------------------------------
    % STEP 4: Sample centers
    %--------------------------------------------------------------------------

    z = cell2mat(z_by_year);
    ix = accumarray(z', 1:length(z), [], @(x){x});
    for i = 1:length(ix)
        if ~isempty(ix{i})
            centers(i,:) = mean(vecs(ix{i}, :));
        end
    end
    
    fprintf('iteration %d done\n', iter)
end
            
            












