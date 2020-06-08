%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script for participation coefficient analysis %%%%%%%%%%%%%%%%%%%%%%%%%
%%% + visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace

clear
clc
close all


%% set paths and directories

addpath(genpath('D:\Mary\work\Toolbox\BCT'));
addpath(genpath('D:\Mary\work\VizStuff'));

dir_comm = 'D:\Mary\work\Lifespan\Data\Communities';
dir_net = 'D:\Mary\work\Lifespan\Data\Network_bootstrap';
dirICN = 'D:\Mary\work\Lifespan\Data';
outdir_fig = 'D:\Mary\work\Lifespan\Figures';


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

ICN = loadname(fullfile(dirICN, 'ICN'));
% ICN: struct with 2 fields --> .labels: cell [15x1] countaining labels of
%                                        the Yeo functional areas
%                               .nodes: vector [114x1] each node has a
%                                       label (from 1 to 15) referring to one
%                                       of the Yeo functional areas

ICN = loadname(fullfile(dirICN, 'ICN'));
% ICN: struct with 2 fields --> .labels: cell [15x1] countaining labels of
%                                        the Yeo functional areas
%                               .nodes: vector [114x1] each node has a
%                                       label (from 1 to 15) referring to one
%                                       of the Yeo functional areas


%% compute participation coefficient for each layer and each iteration

[N, L, ~, ~, iter] = size(comm);

gamma = [0.5, 1, 2];
omega = [0.1, 0.5, 1, 5];

% select omega=0.5 and gamma=2
w = 2;
g = 3;


for l=1:L
    for it=1:iter
        
        pc(:,l,it) = participation_coef(net(:,:,l,it), comm(:,l,w,g,it));
        nullpc(:,l,it) = participation_coef(null_net(:,:,l,it), null_comm(:,l,w,g,it));
        
    end
    clear it
end
clear l


%% visualize: participation coefficient vs age

vec2 = 1:2:39;
xlab2 = 8:2:85;
xlab2(end) = 85;
for i=1:length(vec2)
    lab2{i} = mat2str(xlab2(vec2(i)));
end
clear i

mean_pc = squeeze(mean(pc, 1));
mean_nullpc = squeeze(mean(nullpc, 1));

matrix = cat(3, mean_pc, mean_nullpc);

colors = cbrewer('qual','Paired',12,'pchip');
colLines = colors([2,6],:);
colBounds = colors([1,5],:);
col2plot = cat(3, colLines, colBounds);

[l, b] = plot_BoundedLines(matrix,col2plot);
grid on
legend([b{1}, b{2}],{'Part. Coeff. Observed','Part. Coeff. Null Model'},...
    'Location','best', 'FontSize', 12)
xlabel('Age', 'FontSize', 12)
ylabel('Participation Coefficient', 'FontSize', 12)
xticks(1:2:39)
xticklabels(lab2)


%% correlation between pc and age for each node

var = 1:size(pc,2);

for it=1:iter
    for n=1:N
        
        [rho_pc(n,it), pval_pc(n,it)] = corr(var', pc(n,:,it)', 'type', 'Spearman');
        if pval_pc(n,it) > 0.05
           rho_pc(n,it) = 0; 
        end
        
    end
    clear l
end
clear it

mean_rho = mean(rho_pc, 2);


% visualize on the cortex surface

path = 'D:\Mary\work\VizStuff\plot_on_cortex';
colors = cbrewer('div','PiYG',201,'pchip');
colors = flipud(colors);
col2plot = get_proportionalDivergingColors(mean_rho, colors);

plot_CommOnSurface(path, [1:114]', 'no', col2plot)


%% participation coefficient in the early and late lifespan

% reference to layers 1, 16, 39

var1 = 1:16;
var2 = 1:(39-16);

for it=1:iter
    for n=1:N
        
        [rho_pc116(n,it), pval_pc116(n,it)] = corr(var1', pc(n,1:16,it)',...
            'type', 'Spearman');
        if pval_pc116(n,it) > 0.05
           rho_pc116(n,it) = 0; 
        end
        
        [rho_pc1739(n,it), pval_pc1739(n,it)] = corr(var2', pc(n,17:39,it)',...
            'type', 'Spearman');
        if pval_pc1739(n,it) > 0.05
           rho_pc1739(n,it) = 0; 
        end
        
    end
end

mean_rho116 = mean(rho_pc116,2);
mean_rho1739 = mean(rho_pc1739,2);


% visualization on the cortex surface

colors = cbrewer('div','PiYG',201,'pchip');
colors = flipud(colors);

[col2plot4all, limCmap] = get_proportionalDivergingColors(...
    cat(2, mean_rho, mean_rho116, mean_rho1739), colors);

ageint = {'all', 'young', 'old'};
for idx=1:3
    plot_CommOnSurface(path, [1:114]', 'no', col2plot4all(:,:,idx))
end
clear idx


%% average on Yeo functional areas


clust_icn = ICN.nodes;
lab_icn = ICN.labels;


% sep bars
for i=1:length(lab_icn) 
    nn = find(clust_icn==i);
    
    tmpAll = rho_pc(nn,:);
    tmp_all_pos = tmpAll(find(tmpAll>0)); tmp_all_neg = tmpAll(find(tmpAll<0));
    rho_icn_all_pos(i) = mean(tmp_all_pos(:));
    rho_icn_all_neg(i) = mean(tmp_all_neg(:));
    tmpY = rho_pc116(nn,:);
    tmp_y_pos = tmpY(find(tmpY>0)); tmp_y_neg = tmpY(find(tmpY<0));
    rho_icn_y_pos(i) = mean(tmp_y_pos(:));
    rho_icn_y_neg(i) = mean(tmp_y_neg(:));
    tmpO = rho_pc1739(nn,:);
    tmp_o_pos = tmpO(find(tmpO>0)); tmp_o_neg = tmpO(find(tmpO<0));
    rho_icn_o_pos(i) = mean(tmp_o_pos(:));
    rho_icn_o_neg(i) = mean(tmp_o_neg(:));
    
    rho_icn_all(i) = mean(tmpAll(:));
    rho_icn_y(i) = mean(tmpY(:));
    rho_icn_o(i) = mean(tmpO(:));
end
greens = flipud(colors(1:100,:)); pinks = colors(102:201,:);

max_pos = max([max(rho_icn_all(rho_icn_all>0)),...
    max(rho_icn_y(rho_icn_y>0)),...
    max(rho_icn_o(rho_icn_o>0))]);

max_neg = max([max(abs(rho_icn_all(rho_icn_all<0))),...
    max(abs(rho_icn_y(rho_icn_y<0))),...
    max(abs(rho_icn_o(rho_icn_o<0)))]);

for i=1:length(lab_icn)
    if rho_icn_all(i)>0
        ind_icn_all = round((rho_icn_all(i)*100)/max_pos);
        if isnan(ind_icn_all) || ind_icn_all == 0
            ind_icn_all = 1;
        end
        col_icn_all(i,:) = pinks(ind_icn_all,:);
    else
        ind_icn_all = round((abs(rho_icn_all(i))*100)/max_neg);
        if isnan(ind_icn_all) || ind_icn_all == 0
            ind_icn_all = 1;
        end
        col_icn_all(i,:) = greens(ind_icn_all,:);
    end
    
    if rho_icn_y(i)>0
        ind_icn_y = round((rho_icn_y(i)*100)/max_pos);
        if isnan(ind_icn_y)
            ind_icn_y = 1;
        end
        col_icn_y(i,:) = pinks(ind_icn_y,:);
    else
        ind_icn_y = round((abs(rho_icn_y(i))*100)/max_neg);
        if isnan(ind_icn_y)
            ind_icn_y = 1;
        end
        col_icn_y(i,:) = greens(ind_icn_y,:);
    end
    
    if rho_icn_o(i)>0
        ind_icn_o = round((rho_icn_o(i)*100)/max_pos);
        if isnan(ind_icn_o)
            ind_icn_o = 1;
        end
        col_icn_o(i,:) = pinks(ind_icn_o,:);
    else
        ind_icn_o = round((abs(rho_icn_o(i))*100)/max_neg);
        if isnan(ind_icn_o)
            ind_icn_o = 1;
        end
        col_icn_o(i,:) = greens(ind_icn_o,:);
    end
end

% viualize bar plot: participation coefficient computed on the entire lifespan
% averaged on the Yeo functional areas

figure;
bar_all = bar(rho_icn_all);
bar_all.EdgeColor = 'none';
bar_all.FaceColor = 'flat';
bar_all.CData = col_icn_all;
set(gca, 'XTickLabel', lab_icn, 'XTickLabelRotation', 45, 'FontSize', 12)
ylabel('Correlation: part.coeff. vs age', 'FontSize', 12)
grid on
print(gcf, fullfile(outdir_fig, 'partcoeff_funcareas_all'), '-dtiffn', '-r300')

% viualize bar plot: participation coefficient computed on the early (left bars) 
% and late (right bars) lifespan averaged on the Yeo functional areas

figure;
bar_yo = bar([rho_icn_y', rho_icn_o']);
bar_yo(1,1).EdgeColor = 'none';
bar_yo(1,2).EdgeColor = 'none';
bar_yo(1,1).FaceColor = 'flat';
bar_yo(1,1).CData = [col_icn_y];
bar_yo(1,2).FaceColor = 'flat';
bar_yo(1,2).CData = [col_icn_o];
set(gca, 'XTickLabel', lab_icn, 'XTickLabelRotation', 45, 'FontSize', 12)
ylabel('Correlation: part.coeff. vs age', 'FontSize', 12)
grid on
set(gcf, 'OuterPosition', [481.0000  334.6000  727.2000  509.6000])
print(gcf, fullfile(outdir_fig, 'partcoeff_funcareas_youngold'), '-dtiffn', '-r300')


%  % bar plot referred only to the early lifespan
% figure;
% bar_y = bar(rho_icn_y);
% bar_y.EdgeColor = 'none';
% bar_y.FaceColor = 'flat';
% bar_y.CData = col_icn_y;
% set(gca, 'XTickLabel', lab_icn, 'XTickLabelRotation', 45,...
%     'YLim', [-0.3 0.5], 'FontSize', 12)
% ylabel('Correlation: part.coeff. vs age', 'FontSize', 12)
% grid on


%  % bar plot referred only to the late lifespan
% figure;
% bar_o = bar(rho_icn_o);
% bar_o.EdgeColor = 'none';
% bar_o.FaceColor = 'flat';
% bar_o.CData = col_icn_o;
% set(gca, 'XTickLabel', lab_icn, 'XTickLabelRotation', 45,...
%     'YLim', [-0.3 0.5], 'FontSize', 12)
% ylabel('Correlation: part.coeff. vs age', 'FontSize', 12)
% grid on

