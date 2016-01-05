[~,data] = data_generate();
count = histcounts(data);

g0 = gamrnd(count + 5, 1);
g0 = g0 / sum(g0);

weights = zeros(1, 1000);
for i = 1:1000
    g1 = dpDisrnd(1, g0);
    weights(i) = prod(g1(data(1,:)));
end