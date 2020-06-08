function [h, coord] = myViolinPlot(data,colors)

% INPUT
% data:         n°Obs x n°Var
% colors:       n°Var x 1 cells, each cell contains a 2x3 matrix 

% OUTPUT
% h:        figure properties
% coord:    coord where the density distributions lay

nVar = size(data,2);

base = 1;
coord(1) = base;

% parameter of interest computation
for n=1:nVar
    [f(:,n), y(:,n)] = ksdensity(data(:,n));
    means(n) = mean(data(:,n));
    perc(n,:) = quantile(data(:,n), [0.25 0.75]);
end
clear n
space = max(max(f))/3;

% plot
figure;
for n=1:nVar
    
    % current color
    col = colors{n};
    
    % density plot
    vec1L = cat(1, base-f(:,n), base*ones(length(f(:,n)),1));   % x coord
    vec2L = cat(1, y(:,n), y(:,n));                             % y coord
    h{1,n} = fill(vec1L, vec2L, col(1,:));
    h{1,n}.LineWidth = 1;
    h{1,n}.EdgeColor = col(1,:);
    hold on
    vec1R = cat(1, f(:,n)+base, base*ones(length(f(:,n)),1));   % x coord
    vec2R = cat(1, y(:,n), y(:,n));                             % y coord
    h{2,n} = fill(vec1R, vec2R, col(1,:));
    h{2,n}.LineWidth = 1;
    h{2,n}.EdgeColor = col(1,:);
    hold on
    
    % mean and line plot
    h{3,n} = scatter(base,means(n),50,'MarkerFaceColor',col(2,:),...
        'MarkerEdgeColor',col(2,:));
    hold on
    h{4,n} = line([base base], [perc(n,1) means(n)],...
        'Color', col(2,:), 'LineWidth', 1);
    hold on
    h{5,n} = line([base base], [means(n) perc(n,2)],...
        'Color', col(2,:), 'LineWidth', 1);
    
    if n<nVar
        base = base + max(f(:,n)) + space + max(f(:,n+1));
        coord(n+1) = base;
    end
end
clear n

set(gca,'XTick',coord)
