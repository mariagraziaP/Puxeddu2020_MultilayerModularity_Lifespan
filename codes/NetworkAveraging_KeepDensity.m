function averaged_net = NetworkAveraging_KeepDensity(networks)

% !! NB only undirected networks !!

for d=1:size(networks,3)
    dens(d) = density_und(networks(:,:,d));
end
meanDens = mean(dens);

N = size(networks,1);
nConn = round(meanDens*((N^2-N)/2));

temp_aver_net = mean(networks,3);

masknet = zeros(N,N);
for n=1:size(networks,3)
    bin = networks(:,:,n);
    bin(bin>0) = 1;
    masknet = masknet+bin;
end
clear n

uppermask = triu(masknet);
pos = find(uppermask);
valpos = uppermask(find(uppermask));
[ordval, ordpos] = sort(valpos,'descend');
ordpos = ordpos(1:nConn);
pos = pos(ordpos);

averaged_net = zeros(N,N);
averaged_net(pos) = temp_aver_net(pos);
averaged_net = averaged_net + averaged_net';