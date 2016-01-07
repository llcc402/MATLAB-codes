clear
clc

[distro, distro_noise, data, kl] = data_generate();

tic;
[distro_smooth, G0] = synthc_particle(data, 1, 2000, 100, 1, 1, 0);
toc