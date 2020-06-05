%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script to build the ensemble of equivalent multilayer networks %%%%%%%%
%%% and visualize networks distribution across years with the binning %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace

clear
clc
close all

%% set paths and directories

dir_data = 'D:\Mary\work\Lifespan\Data';
savedir_net = 'D:\Mary\work\Lifespan\Data\Network_bootstrap';


%% load data

load(fullfile(dir_data,'MLNetwork')) %network
% load a matrix of dimension [N*N*T]: N=114(nodes), L=620(subjects)
% [after the preprocessing we concatenated all the anatomical networks (one
% for each subject) in a 3D matrix from the youngest subject to the oldest] 

loadname(fullfile(dir_data,'Age')) %age
% vector [L*1] containing ages associated to each network


%% binning of the dataset + multilayer networks construction

[N, ~, T] = size(network);

binwidth = 2;       % range of years for each bin
iter = 1000;        % resampling iterations
sampsize = 10;      % number of networks sampled for each bin

bound = min(age):binwidth:(max(age)-1);

net = zeros(N, N, length(bound), iter);

for nb=1:length(bound)
    
    % binning
    age1 = bound(nb);
    if nb==length(bound)
        age2 = bound(nb+2);
    else
        age2 = bound(nb+1);
    end
    
    lower_bound = find(age==age1);
    higher_bound = find(age==age2);
    
    idx = lower_bound(1):higher_bound(end);
    
    % resampling
    for it=1:iter
        
        curr_idx = idx(randi([1 length(idx)], 1, 10));
        curr_net = network(:,:,curr_idx);
        net(:,:,nb,it) = NetworkAveraging_KeepDensity(curr_net);
        
    end
    
end

save(fullfile(savedir_net, 'Network_bin2years'), 'net')


%% visualize histogram with age distribution + binning 

NL = length(bound);     % number of layers
colors = cbrewer('qual', 'Set3', NL, 'pchip');

ind = randperm(NL);
colors = colors(ind, :);

figure;

for l=1:NL
    
    vec = [];
    
    if l==NL
        num = binwidth+1;
    else
        num = binwidth;
    end
    
    for idx=1:num
        age1 = bound(l)+idx-1;
        vec = [vec; repmat(age1, [21,1])];
        clear age1
    end
        
    hh = histogram(vec,...
        'BinLimits', [0,90],...
        'BinWidth', 1,...
        'EdgeColor', 'none',...
        'FaceColor', colors(l,:));
    hold on
    
    clear vec idx num
end

hold on

histogram(age,...
    'BinLimits', [7,86],...
    'BinWidth', 1,...
    'EdgeColor', 'k',...
    'FaceColor', [0, 0.4470, 0.7410],...
    'FaceAlpha', 1)

axis equal
axis tight

ylim([0,21])
xlim([7, 86])

set(gca, 'FontSize', 12)
xlabel('Age', 'FontSize', 12)

title('bin = 2years, NG = 39 (num. layers)', 'FontSize', 12)

hold off
