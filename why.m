tic;
K = 1000;
kl_s = zeros(1, K);
sum_g2 = zeros(1,K);
G0_hard = [];
G1_hard = [];
G2_hard = [];
parfor k = 1:K
    G0 = gem(100, 5);
    G1 = dpDisrnd(1, G0);
    
    n = 0;
    G2 = smoothSample(G0,G1);
    kl_s(k) = symKL(G1, G2); 
    if kl_s(k) > 10
        G0_hard = [G0_hard; G0];
        G1_hard = [G1_hard; G1];
        G2_hard = [G2_hard; G2];
    end
    
    sum_g2(k) = sum(G2);
    fprintf(['iter ', num2str(k),' done \n'])
end
toc
