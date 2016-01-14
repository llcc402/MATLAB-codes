clear
clc

%--------------------------------------------------------------------------
% STEP 1: test kl with different bounds
%--------------------------------------------------------------------------
% kl = zeros(1000,11);
% activated = 100;
% gamma = 5; % the concentration of G0
% alpha = 1; % the concentration of G(1:end)
% tic;
% for i = 1:1000
%     for j = 1:10
%         G0 = gem(activated, gamma);
%         G1 = dpDisrnd(alpha, G0);
%         G2 = smoothSample(G0, G1, j, alpha);
%         kl(i,j) = symKL(G1, G2);
%     end
% end
% fprintf('The time of smoothed sampling for 10000 times is: %f\n', toc)
% 
% tic
% for i = 1:1000
%     G0 = gem(activated, gamma);
%     G1 = dpDisrnd(alpha, G0);
%     G2 = dpDisrnd(alpha, G0);
%     kl(i, 11) = symKL(G1, G2);
% end
% fprintf('The time of direct sampling for 1000 times is: %f\n', toc)
% boxplot(kl)
% title('The symetric kl-divergence with different bounds')


%--------------------------------------------------------------------------
% STEP 2: post sample with synthetic data
%--------------------------------------------------------------------------
clear
clc

activateN = 100;
particleN = 2000;
gibsN = 10;
alpha = 1;
B = 1;
gamma = 0.001;
[distro, distro_noise, data, kl] = data_generate();

% smooth sampling
tic; 
[distro_smooth, G0] = synthc_particle(data, gibsN, particleN,...
    activateN, alpha, B, gamma);
fprintf('The time of smooth sampling for %d iterations is: %f seconds\n', gibsN, toc)

% traditional sampling
tic;
distro_trad = traditional_particle(data, particleN, activateN, gamma,...
    alpha, gibsN);
fprintf('The time of traditional sampling for %d iterations is: %f seconds\n', gibsN, toc)

