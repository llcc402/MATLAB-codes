% This file is used to compare the kl-divergence between the traditional dp
% sampling method with our smoothed dp sampling method.

tic;
N = 1000;
B = 3;
alpha = 1;
gamma = 5;

kl_tradition = zeros(1, N);
kl_smooth = zeros(1, N);

for i = 1:N
    
    G0 = gem(100, gamma);
    G1 = dpDisrnd(alpha, G0);
    
    G2 = dpDisrnd(alpha, G0);
    kl_tradition(i) = symKL(G1, G2);
    
    G2 = smoothSample(G0, G1, B, 1);
    kl_smooth(i) = symKL(G1, G2);
    
end
kl_both = [kl_smooth', kl_tradition'];
boxplot(kl_both)

toc
