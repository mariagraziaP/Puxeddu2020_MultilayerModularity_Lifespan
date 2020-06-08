%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script to compute the single-layer modularity %%%%%%%%%%%%%%%%%%%%%%%%%
%%% + visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace

clear
clc
close all


%% set paths and directories

dir_net = 'D:\Mary\work\Lifespan\Data\Network_bootstrap';
dir_comm = 'D:\Mary\work\Lifespan\Data\Communities';
savedir = 'D:\Mary\work\Lifespan\Data\Indices';
outdir_fig = 'D:\Mary\work\Lifespan\Figures';


%% load data

load(fullfile(dir_net,'Network_bin2years')) %net
% matrix of dimension [N*N*L*ITER]:  N=114(nodes), 
%                                    L=39(layers),
%                                    ITER=1000(number of multilayer networks)

load(fullfile(dir_net,'Communities_bin2years')) %comm
% matrix of dimension [N*L*4*3*ITER]:   N=114(nodes), 
%                                       L=39(layers),
%                                       4= number of omega explored
%                                       3= number of gamma explored
%                                       ITER=1000(number of multilayer networks)

%% compute single-layer modularity

gamma = [0.5, 1, 2];
omega = [0.1, 0.5, 1, 5];

[N, L, ~, ~, iter] = size(comm);

slmod = zeros(L, iter, length(omega), length(gamma));
slmod_null = slmod;

for w=1:length(omega)
    for g=1:length(gamma)
        for s=1:L
            for it=1:iter
                
                slmod(s,it,w,g) = Modularity(...
                    net(:,:,l,it),...
                    comm(:,l,w,g,it),...
                    gamma(g), 'und');
                
                slmod_null(s,it,w,g) = Modularity(...
                    null_net(:,:,l,it),...
                    null_comm(:,l,w,g,it),...
                    gamma(g), 'und');
            end
        end
    end
end

save(fullfile(savedir_cn, 'SLModularity_bin2years'), 'slmod')
save(fullfile(savedir_cn, 'SLModularity_NullMod_bin2years'), 'slmod_null')


%% visualization only for omega=0.5 and gamma=2

w = 2;
g = 3;

%% visualize: single-layer modularity vs age

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

matrix = cat(3, squeeze(slmod(:,:,w,g)), squeeze(slmod_null(:,:,w,g)));

[l, b] = plot_BoundedLines(matrix, col2plot);
grid on
legend([b{1}, b{2}],{'Q observed', 'Q null model'},...
    'Location', 'best', 'FontSize', 12)
xlabel('Age', 'FontSize', 12)
xticks(1:2:39)
xticklabels(lab2)


%% number/weights of intra-cluster links vs age

conn_intra = zeros(L, iter);
weights_intra = conn_intra;

for s=1:L
    for it=1:iter
        
        curr_net = squeeze(net(:,:,s,it));
        curr_comm = squeeze(comm(:,s,w,g,it));
        
        tot_conn = length(find(curr_net));
        lab = unique(curr_comm);
        
        count = zeros(length(lab));
        weights = count;
        for jdx=1:length(lab)
            indcl = find(curr_comm==lab(jdx));
            count(jdx) = length(find(curr_net(indcl, indcl)));
            weights(jdx) = sum(sum(curr_net(indcl, indcl)));
        end
        
        conn_intra(s, it) = sum(count);
        weights_intra(s, it) = sum(weights);
        
    end
end


% visualization
colors = cbrewer('qual','Paired',12,'pchip'); 
colLines = colors([4, 10],:);
colBounds = colors([3, 9],:);
col2plot = cat(3, colLines, colBounds);

matrix = cat(3, conn_intra, weights_intra);

fig = figure;

left_color = colors(4,:);
right_color = colors(10,:);
set(fig, 'defaultAxesColorOrder', [left_color; right_color]);
yyaxis left

[l, b] = plot_BoundedLines(conn_intra, col2plot(1,:,:));
grid on
hold on
plot(mean(conn_intra,2), 'Color', colors(4,:), 'LineStyle', '-')

ylabel('Conn. Intra Cluster', 'FontSize', 12)
xlabel('Age', 'FontSize', 12)
xticks(1:2:39)
xticklabels(lab2)
yyaxis right

[l, b] = plot_BoundedLines(weights_intra, col2plot(2,:,:));
ylabel('Weights Intra Cluster', 'FontSize', 12)
hold on
plot(mean(weights_intra,2), 'Color', colors(10,:), 'LineStyle', '-')



% scatterplots
colors = cbrewer('div', 'RdBu', 39, 'pchip');

% modularity vs intra-cluster weights
figure;
hold on
for it=1:size(mod,2)
    scatter(mod(:,it), weights_intra(:,it), 36, colors, 'filled', 'MarkerFaceAlpha', 0.7)
    hold on
end
box on
xlabel('Single Layer Modularity', 'FontSize', 12)
ylabel('Intra Cluster Weights', 'FontSize', 12)


% modularity vs number of intra-cluster connections
figure;
for it=1:size(mod,2)
    scatter(mod(:,it), conn_intra(:,it), 36, colors, 'filled', 'MarkerFaceAlpha', 0.7)
    hold on
end
box on
xlabel('Single Layer Modularity', 'FontSize', 12)
ylabel('Intra Cluster Connections', 'FontSize', 12)

