%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script to compute the clusters number and size for each layer %%%%%%%%%
%%% + visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace

clear
clc
close all


%% set paths and directories

dir_comm = 'D:\Mary\work\Lifespan\Data\Communities';
savedir_cn = 'D:\Mary\work\Lifespan\Data\Indices';
outdir_fig = 'D:\Mary\work\Lifespan\Figures';


%% load data

load(fullfile(dir_net,'Communities_bin2years')) %comm
% matrix of dimension [N*L*4*3*ITER]:   N=114(nodes), 
%                                       L=39(layers),
%                                       4= number of omega explored
%                                       3= number of gamma explored
%                                       ITER=1000(number of multilayer networks)

%% compute clusters number

gamma = [0.5, 1, 2];
omega = [0.1, 0.5, 1, 5];

[N, L, ~, ~, iter] = size(comm);

clust_num = zeros(L, iter, length(omega), length(gamma));
clust_size = zeros(L, iter, length(omega), length(gamma));

for it=1:iter
    for w=1:length(omega)
        for g=1:length(gamma)
            for s=1:L
                curr_comm = comm(:,nl,w,g,it);
                curr_comm_lab = unique(curr_comm);
                
                clust_num(s,it,w,g) = length(curr_comm_lab);
                clust_size_all = zeros(length(curr_comm_lab), 1);
                for idx=1:length(curr_comm_lab)
                    clust_size_all(idx) = length(find(curr_comm==curr_comm_lab(idx)));
                end
                clust_size(s,it,w,g) = mean(clust_size_all);
            end
        end
    end
end

save(fullfile(savedir_cn, 'ClusterNumber_bin2years'), 'clust_num')
save(fullfile(savedir_cn, 'ClusterSize_bin2years'), 'clust_size')


%% visualize: violin plots

% only for omega=0.5

% average across layers
cn_allL = squeeze(mean(clust_num, 1));
cs_allL = squeeze(mean(clust_size, 1));

% select omega .5
cn2plot = squeeze(cn_allL(:,2,:)); % 1000x3 (iter x gammas)
cs2plot = squeeze(cs_allL(:,2,:));

% get colors
[cb] = cbrewer('qual','Paired',12,'pchip');


% make violins: cluster number
colcn{1} = cb([2 1],:); colcn{2} = cb([2 1],:); colcn{3} = cb([2 1],:);

[h, coord] = myViolinPlot(cn2plot, colcn);

grid on
set(gca, 'XTickLabel',{'0.5','1','2'}, 'FontSize', 12)
xlabel('\gamma', 'FontSize', 16)
ylabel('Average Clusters Number', 'FontSize', 14)
set(gca, 'XLim', [-1 9])

print(gcf, fullfile(outdir_fig, 'violinClustNum_bin2year.tif'), '-dtiffn', '-r300')

% make violins: cluster size
colcs{1} = cb([6 5],:); colcs{2} = cb([6 5],:); colcs{3} = cb([6 5],:);

[h, coord] = myViolinPlot(cs2plot, colcs);

grid on
set(gca, 'XTickLabel',{'0.5','1','2'}, 'FontSize', 12)
xlabel('\gamma', 'FontSize', 16)
ylabel('Average Clusters Size', 'FontSize', 14)
set(gca, 'YLim', [0 62])
xlim([0.5 5])

print(gcf, fullfile(outdir_fig, 'violinClustSize_bin2year.tif'), '-dtiffn', '-r300')



%% visualize: scatterplot

col1 = cbrewer('qual', 'Set1', 6, 'pchip');
col2 = cbrewer('qual', 'Set2', 6, 'pchip');

cback = [col1(3,:); col2(6,:); col1(4,:)];

colors = cbrewer('div', 'RdBu', 39, 'pchip');

back = repmat([15.5, 19.5, 25], 15,1);

for w=2 %:size(clust_num,3)
    
    figure;
    a = area(back);
    a(1).FaceColor = cback(3,:);
    a(2).FaceColor = cback(2,:);
    a(3).FaceColor = cback(1,:);
    a(1).FaceAlpha = 0.6;
    a(2).FaceAlpha = 0.6;
    a(3).FaceAlpha = 0.6;
    a(1).LineStyle = 'none';
    a(2).LineStyle = 'none';
    a(3).LineStyle = 'none';
    hold on
    
    for g=1:size(clust_num,4)
        cn = squeeze(clust_num(:,:,w,g));
        cs = squeeze(clust_size(:,:,w,g));
             
        allx = [];
        ally = [];
        for l=1:layers
            allx = cat(1, allx, cn(l,:)');
            ally = cat(1, ally, cs(l,:)');
        end
        clear l
        uncoord = unique([allx ally], 'rows');
        for n=1:size(uncoord,1)
            [x(:,n), y(:,n)] = RandomCoordCircle(...
                uncoord(n,1), uncoord(n,2), 0.7, layers);
        end
        clear n
        
        for l=1:layers
            tmpcn = cn(l,:)';
            tmpcs = cs(l,:)';
            for n=1:size(uncoord,1)
                count = length(find(tmpcn == uncoord(n,1)));
                sizedata(n,l) = (count*200)/1000;
            end
            clear n
        end
        clear l
        for n=1:size(uncoord,1)
            curr_size = sizedata(n,:);
            [ordsize, ord] = sort(curr_size, 'ascend');
            pos = find(ordsize>0);
            age2plot = ord(pos);
            for q=1:length(age2plot)
                scatter(x(q,n), y(q,n), 'SizeData', ordsize(pos(q)),...
                        'MarkerFaceColor', colors(age2plot(q),:),...
                        'MarkerEdgeColor', colors(age2plot(q),:))
                hold on
            end
            clear l curr_size ordsize ord pos age2plot
        end
        clear n


        clear cn cs allx ally uncoord x y
    end
    grid on
    xlabel('Clusters Number')
    ylabel('Clusters Size')
    xlim([1 15])
    xticks([1:2:15])
end
clear w

print(gcf, fullfile(outDir, 'ClustNum_vs_ClustSize.tif'), '-dtiffn', '-r300')
