% This file is used to compare the kl-divergence between the traditional dp
% sampling method with our smoothed dp sampling method.

tic;
N = 1000;

kl_tradition = zeros(1, N);

parfor i = 1:N
    
    G0 = gem(100, 5);
    G1 = dpDisrnd(1, G0);
    
    G2 = dpDisrnd(1, G0);
    kl_tradition(i) = symKL(G1, G2);
    
end

toc
