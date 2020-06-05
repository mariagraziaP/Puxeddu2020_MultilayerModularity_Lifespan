%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script to run the multilayer modularity optimization %%%%%%%%%%%%%%%%%%
%%% on the ensemble of multilayer networks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace

clear
clc
close all

%% set paths and directories

dir_net = 'D:\Mary\work\Lifespan\Data\Network_bootstrap';
savedir_comm = 'D:\Mary\work\Lifespan\Data\Communities';

%% load data

load(fullfile(dir_net,'Network_bin2years')) %net
% matrix of dimension [N*N*L*ITER]:  N=114(nodes), 
%                                    L=39(layers),
%                                    ITER=1000(number of multilayer networks)

%% 

[N, ~, L, iter] = size(net);

% set genlouvain parameters
gamma = [0.5, 1, 2];        % spatial resolution parameters 
omega = [0.1, 0.5, 1, 5];   % inter-layer coupling
limit = 10000;
verbose = false;

% run modularity optimization on all the multilayer networks
for it=1:iter
    
    Bnet = spalloc(N*L, N*L, N*N*L+2*N*L);
    Bg = Bnet;
    Bw = Bnet;
    
    for s=1:L
        idx = (1:N)+(s-1)*N;
        Bnet(idx, idx) = net(:,:,s,it);
        k = sum(net(:,:,s,it));
        Bg(idx, idx) = (k.'*k)/sum(k);
    end
    Bw((N+1):L*N, 1:N*(L-1)) = eye(N*(L-1));
    Bw = Bw + Bw';
    
    for w=1:length(omega)
        for g=1:length(gamma)
            
            Bfull = (Bnet - gamma(g)*Bg) + omega(w)*Bw;
            S = genlouvain(Bfull, limit, verbose);
            comm(:,:,w,g,it) = reshape(S, [N, L]);
            
        end
    end
end

save(fullfile(savedir_comm, 'Communities_bin2yeasrs'))
