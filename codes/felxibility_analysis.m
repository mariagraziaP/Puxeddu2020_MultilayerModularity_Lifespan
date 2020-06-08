%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script for flexibility analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% + visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace

clear
clc
close all


%% set paths and directories

dir_comm = 'D:\Mary\work\Lifespan\Data\Communities';
savedir = 'D:\Mary\work\Lifespan\Data\Indices';
dirICN = 'D:\Mary\work\Lifespan\Data';


%% load data

load(fullfile(dir_comm,'Communities_bin2years')) %comm
% matrix of dimension [N*L*4*3*ITER]:   N=114(nodes), 
%                                       L=39(layers),
%                                       4= number of omega explored
%                                       3= number of gamma explored
%                                       ITER=1000(number of multilayer networks)
load(fullfile(dir_comm,'Communities_NullModel_bin2years')) %null_comm

ICN = loadname(fullfile(dirICN, 'ICN'));
% ICN: struct with 2 fields --> .labels: cell [15x1] countaining labels of
%                                        the Yeo functional areas
%                               .nodes: vector [114x1] each node has a
%                                       label (from 1 to 15) referring to one
%                                       of the Yeo functional areas


%% compute flexibility

[N, L, ~, ~, iter] = size(comm);

gamma = [0.5, 1, 2];
omega = [0.1, 0.5, 1, 5];

% select omega=0.5 and gamma=2
w = 2;
g = 3;

% flexibility
for it=1:iter
    
    [Fmat(:,:,it), nodF(:,it), layF(it,:)] = myFlexibility(comm(:,:,w,g,it));
    [nullFmat(:,:,it), nullnodF(:,it), nulllayF(it,:)] = myFlexibility(null_comm(:,:,w,g,it));

end


%% visualize layer flexibility vs age

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


matrix = cat(3, layF', nulllayF');


figure;
[l, b] = plot_BoundedLines(matrix, col2plot);
grid on
legend([b{1}, b{2}],{'F observed', 'F null model'},...
    'Location', 'best',  'FontSize', 12)
xlabel('Age', 'FontSize', 12)
xticks(1:2:39)
xticklabels(lab2)
ylabel('Layer Flexibility', 'FontSize', 12)


%% compute correlation laver flexibility vs age

mean_lf = mean(layF,1);
var = 1:size(matrix,1);
[rho_fl, pval_fl] = corr(var', mean_lf', 'type', 'Spearman');



%% visualize nodes flexibility

colors = flipud(cbrewer('seq', 'Reds', 100, 'pchip'));

mean_nodF = mean(nodF, 2);

path = 'D:\Mary\work\VizStuff\plot_on_cortex';

maxNF = max(mean_nodF);
minNF = min(mean_nodF);
diffNF = maxNF-minNF;
col2plot = zeros(size(nodF,1),3);
for n=1:size(nodF,1)
    ind = round(((mean_nodF(n)-minNF)*size(colors,1))/diffNF); 
    if ind==0
        ind = 1;
    end
    col2plot(n,:) = colors(ind,:);
end
clear n

plot_CommOnSurface(path, [1:114]', 'no', col2plot)


%% node flexibility in the early and late lifespan

% reference to layers 1, 16, 39

for it=1:iter
    
    [~, nodF116(:,it), ~] = myFlexibility(comm(:,1:16,w,g,it));
    [~, nodF1739(:,it), ~] = myFlexibility(comm(:,17:end,w,g,it));

end

mean_nodF116 = mean(nodF116,2);
mean_nodF1739 = mean(nodF1739,2);

% plot

ageint = {'all', 'young', 'old'};

[col2plot4all, limCmap] = get_proportionalColors(...
    cat(2, mean_nodF, mean_nodF116, mean_nodF1739),...
    colors);


for idx=1:size(col2plot4all,3)
    
    plot_CommOnSurface(path, [1:114]', 'no', col2plot4all(:,:,idx))
    print(gcf, fullfile(outDir, sprintf('NodeFlex_%s', ageint{idx})), '-dtiffn', '-r300')

end
clear idx


%% average on Yeo functional areas


clust_icn = ICN.nodes;
lab_icn = ICN.labels;

t = 1;
for i=1:length(lab_icn) 
    
    nn = find(clust_icn==i);
    tmpAll = nodF(nn,:);
    tmpY = nodF116(nn,:);
    tmpO = nodF1739(nn,:);
    
    flex_icn_all(i) = mean(tmpAll(:));
    flex_icn_y(i) = mean(tmpY(:));
    flex_icn_o(i) = mean(tmpO(:));
    
    flex_icn_yo(t) = flex_icn_y(i);
    flex_icn_yo(t+1) = flex_icn_o(i);
    t = t+2;
end

% viualize bar plot: node flexibility computed on the entire lifespan
% averaged on the Yeo functional areas

figure;
for i=1:length(lab_icn)
    col_icn_all = round((flex_icn_all(i)*1000)/max([max(flex_icn_y), max(flex_icn_o), max(flex_icn_all)]));
    
    rectangle('Position', [i, 0, 0.8, flex_icn_all(i)],...
        'FaceColor', colors(col_icn_all,:), 'EdgeColor', 'none')
    
    hold on
end
grid on
box on
set(gca, 'XTick', unique(clust_icn)+0.2,...
    'XTickLabel', lab_icn,...
    'XTickLabelRotation', 45,...
    'XLim', [0.8 16],...
    'FontSize', 12)
ylabel('Node Flexibility', 'FontSize', 12)


% viualize bar plot: node flexibility computed on the early (left bars) and late
% (right bars) lifespan averaged on the Yeo functional areas

figure;
c = 1;
b = 1;
for i=1:2*length(lab_icn)
    col_icn_o = round((flex_icn_yo(i)*1000)/max([max(flex_icn_y), max(flex_icn_o), max(flex_icn_all)]));
    
    rectangle('Position', [b, 0, 0.8, flex_icn_yo(i)],...
        'FaceColor', colors(col_icn_o,:), 'EdgeColor', 'none')
    clear base
    hold on
    if mod(i,2)==0
        b = b+2;
    else
        b = b+1;
        xt(c) = b;
        c = c+1;
    end
end
grid on
box on
set(gca, 'XTick', xt(1:15),...
    'XTickLabel', lab_icn,...
    'XTickLabelRotation', 45,...
    'XLim', [0.8 45],...
    'FontSize', 12)
ylabel('Node Flexibility', 'FontSize', 12)
