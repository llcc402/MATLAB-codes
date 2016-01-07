% This file is used to test how the kl-divergece change with the change of
% the bound from 1 to 10 with pace 1.
clear
clc

N = 100;
K = 10;
kl_bound_change = zeros(N, K);
        
tic;
parfor k = 1:K
    temp = zeros(N,1);
    for i = 1:N
        G0 = gem(100, 5);
        G1 = dpDisrnd(1, G0);
        
        G2 = smoothSample(G0, G1, k);
        temp(i) = symKL(G1, G2);
    end
    kl_bound_change(:, k) = temp;
end
toc

boxplot(kl_bound_change)
title('The kl-divergence with the bound change from 1 to 10.')