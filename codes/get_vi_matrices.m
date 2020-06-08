%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script to compute the variation of information %%%%%%%%%%%%%%%%%%%%%%%%
%%% + visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace

clear
clc
close all


%% set paths and directories

dir_comm = 'D:\Mary\work\Lifespan\Data\Communities';
savedir_vi = 'D:\Mary\work\Lifespan\Data\Indices';


%% load data

load(fullfile(dir_net,'Communities_bin2yeasrs')) %comm
% matrix of dimension [N*L*4*3*ITER]:   N=114(nodes), 
%                                       L=39(layers),
%                                       4= number of omega explored
%                                       3= number of gamma explored
%                                       ITER=1000(number of multilayer networks)

%%

gamma = [0.5, 1, 2];
omega = [0.1, 0.5, 1, 5];

[N, L, ~, ~, iter] = size(comm);

aux = ones(L, L);

[ind1,ind2] = ind2sub(size(aux), find(triu(aux,1)==1));
ncomp = length(ind1);

% compute variation of information
vi = zeros(L, L, length(omega), length(gamma), iter);
for it=1:iter
    for w=1:length(omega)
        for g=1:length(gamma)
            temp = zeros(L, L);
            for t=1:ncomp
                temp(ind1(t),ind2(t)) = ...
                    partition_distance(...
                    comm(:, ind1(t), w, g, it),...
                    comm(:, ind2(t), w, g, it));
            end
            vi(:,:,w,g,it) = temp + temp';
        end
    end
end

save(fullfile(savedir_vi ,'VI_bin2years'),'vi')



%% visualization of the VI matrices for each gamma and omega

mean_vi = mean(vi, 5);

vec2 = 0:10:39;
vec2(1) = 1;
vec2(end+1) = 39;
xlab2 = 8:2:85;
xlab2(end) = 85;
for i=1:length(vec2)
    lab2{i} = mat2str(xlab2(vec2(i)));
end
clear i

map = parula(400);

for w=1:length(omega)
    for g=1:length(gamma)
        
        figure;
        imagesc(mean_vi(:,:,w,g))
        axis square
        
        xlabel('Age', 'FontSize', 14)
        ylabel('Age', 'FontSize', 14)
        
        set(gca, 'XTick', vec2,...
            'XTickLabel', lab2,...
            'YTick', vec2,...
            'YTickLabel', lab2,...
            'FontSize', 14)

        cl(:,w,g) = get(gca, 'CLim');
        set(gca, 'CLim', [0 0.43])
        
        if w==1
            ylabel({cat(2,...
                '\bf \gamma',...
                sprintf(' = %g', gamma(g)));...
                '\rm Age'},...
                'FontSize', 14)
        end
        if g==1
            title(cat(2,...
                '\omega',...
                sprintf(' = %g', omega(w))),...
                'FontSize', 14)
        end
        
        print(gcf, fullfile(savedir_vi,...
            sprintf('VImatrix_bin2year_w%gg%g.tif', omega(w), gamma(g))),...
            '-dtiffn', '-r300')
    
    end
end

    
