%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script to observe how the topographic distribution of modules %%%%%%%%%
%%%%%%%%%%%% on the cortex changes across the lifespan %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% + visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace

clear
clc
close all


%% set paths and directories

dir_comm = 'D:\Mary\work\Lifespan\Data\Communities';
dir_net = 'D:\Mary\work\Lifespan\Data\Network_bootstrap';


%% load data

load(fullfile(dir_comm,'Communities_bin2years')) %comm
% matrix of dimension [N*L*4*3*ITER]:   N=114(nodes), 
%                                       L=39(layers),
%                                       4= number of omega explored
%                                       3= number of gamma explored
%                                       ITER=1000(number of multilayer networks)
load(fullfile(dir_comm,'Communities_NullModel_bin2years')) %null_comm

load(fullfile(dir_net,'Network_bin2years')) %net
% matrix of dimension [N*N*L*ITER]:  N=114(nodes), 
%                                    L=39(layers),
%                                    ITER=1000(number of multilayer networks)
load(fullfile(dir_net,'NullNetwork_bin2years')) %null_net


%%

[N, L, ~, ~, iter] = size(comm);

gamma = [0.5, 1, 2];
omega = [0.1, 0.5, 1, 5];

% select omega=0.5 and gamma=2
w = 2;
g = 3;

% nodes from left and right hemispheres
lh = 1:57;
rh = 58:114;


% compute distribution on hemispheres across lifespan
ncl_comm = zeros(L, iter);
ncl_null_comm = ncl_comm;

for l=1:L
    for it=1:iter
        
        curr_ris = squeeze(comm(:, l, w, g, it));
        ncl_comm(l,it) = ...
            length(find(...
            ismember(unique(curr_ris(lh)), unique(curr_ris(rh)))));
        
        curr_nullris = squeeze(null_comm(:, l, w, g, it));
        null_ncl_comm(l,it) = ...
            length(find(ismember(...
            unique(curr_nullris(lh)), unique(curr_nullris(rh)))));
        
    end
    clear it
end
clear l


%% visualize: number of modules spanning the two hemispheres vs age

vec2 = 1:2:39;
xlab2 = 8:2:85;
xlab2(end) = 85;
for i=1:length(vec2)
    lab2{i} = mat2str(xlab2(vec2(i)));
end
clear i

colors = cbrewer('qual','Paired',12,'pchip');
colLines = colors([2, 6],:);
colBounds = colors([1, 5],:);
col2plot = cat(3, colLines, colBounds);

matrix = cat(3, ncl_comm, null_ncl_comm);

[l, b] = plot_BoundedLines(matrix, col2plot);
grid on
legend([b{1}, b{2}],{'Observed','Null Model'},...
    'Location','best', 'FontSize', 12)
xlabel('Age', 'FontSize', 12)
xticks(1:2:39)
xticklabels(lab2)
ylabel('Num of Clusters spanning both Hemispheres', 'FontSize', 12)


%% find which nodes belong to the last module

nodes = [];

for it=1:iter
    
    curr_comm = comm(:, 39, w, g, it);
    
    lablh = unique(curr_comm(lh));
    labrh = unique(curr_comm(rh));
    comm_lab = lablh(ismember(lablh, labrh));
    tmp = [];
    for idx=1:length(comm_lab)
        tmp = [tmp; find(curr_comm==comm_lab(idx))];
    end
    nodes = [nodes; tmp];
    
    clear tmp curr_ris idx lablh labrh comm_lab
end
clear it

for n=1:N
    freq(n,1) = length(find(nodes==n));
end
clear n


% visualize through a hit-map on the cortex

path = 'D:\Mary\work\VizStuff\plot_on_cortex';

colors = cbrewer('seq', 'Blues', 1000, 'pchip');
[col2plot limCmap] = get_proportionalColors(freq, colors);

plot_CommOnSurface(path, [1:114]', 'no', col2plot)


%% relationship between topographic distribution and inter-hemispheric links

%% compute the number and weights of the intra/inter-hemispehric links across the lifespan

for l=1:L
    for it=1:iter
        
        curr_net = net(:,:,l,it);
        
        dens(l,it) = density_und(net);
        
        wei_lh(l,it) = sum(sum(net(lh,lh)));
        wei_rh(l,it) = sum(sum(net(rh,rh)));
        wei_ih(l,it) = 2*sum(sum(net(lh,rh)));
        
        conn_lh(l,it) = length(find(net(lh,lh)));
        conn_rh(l,it) = length(find(net(rh,rh)));
        conn_ih(l,it) = 2*length(find(net(lh,rh)));
        
    end
end

% compute correlations

var = 1:size(conn_ih,1);

[rho_wam_s, pval_wam_s] = corr(var', mean(wei_ih,2), 'type', 'Spearman');
[rho_cam_s, pval_cam_s] = corr(var', mean(conn_ih,2), 'type', 'Spearman');

[rho_wcm_s, pval_wcm_s] = corr(mean(ncl_comm,2), mean(wei_ih,2), 'type', 'Spearman');
[rho_ccm_s, pval_ccm_s] = corr(mean(ncl_comm,2), mean(conn_ih,2), 'type', 'Spearman');


%% visualize bounded lines

colors = cbrewer('qual','Paired',12,'pchip'); 
colLines = colors([4, 10, 8],:);
colBounds = colors([3, 9, 7],:);
col2plot = cat(3, colLines, colBounds);

% number of intra/inter-hemispheric links across lifespanù

matrix_conn = cat(3, conn_lh, conn_rh, conn_ih);

figure;
[l, b] = plot_BoundedLines(matrix_conn, col2plot);
grid on
hold on
ylabel('Number of edges', 'FontSize', 12)
xlabel('Age', 'FontSize', 12)
xticks(1:2:39)
xticklabels(lab2)
legend([b{1}, b{2}, b{3}], {'Left Hem','Right Hem', 'Inter-Hem'},...
    'Location','best', 'FontSize', 12)


% weights of intra/inter-hemispheric links across lifespanù

matrix_wei = cat(3, wei_lh, wei_rh, wei_ih);

figure;
[l, b] = plot_BoundedLines(matrix_wei, col2plot);
grid on
hold on
ylabel('Weights of edges', 'FontSize', 12)
xlabel('Age', 'FontSize', 12)
xticks(1:2:39)
xticklabels(lab2)
legend([b{1}, b{2}, b{3}], {'Left Hem','Right Hem', 'Inter-Hem'},...
    'Location','best', 'FontSize', 12)


%% rain-cloud plot: inter-hemispheric weights vs number of cluster spanning both hemispehres wrt age

cllab = unique(ncl_comm);

colors = cbrewer('div', 'RdBu', 39, 'pchip');

colrain{1} = [0.8 0.8 0.8];
colrain{2} = [0.3 0.3 0.3];

for i=1:length(cllab)
    
    pos = find(ncl_comm==cllab(i));
    wei2vio{i} = wei_ih(pos);
    [age, iter] = ind2sub([39, 1000], pos);
    
    for t=1:39
        pt = find(age==40-t);
        colrain{3}(t,:) = colors(40-t,:);
        if not(isempty(pt))
            weiallage{t,i} = wei_ih(pos(pt));
        else
            weiallage{t,i} = [];
        end
    end
    
end

[h, coord] = myRainCloudPlot_multicolor_ls(wei2vio, weiallage, colrain);
box on
grid on
set(gca, 'XTick', coord,...
    'XTickLabel', {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'},...
    'XLim', [0.9 1.9],...
    'FontSize', 14)

xlabel('Number of Clusters spanning both Hemispheres', 'FontSize', 20)
ylabel('Inter-Hemispheric Weights', 'FontSize', 20)

set(gcf, 'OuterPosition', [481 334.6 888 509.6])
