function [h, coord] = myRainCloudPlot_multicolor_ls(data, subvar, color)

% INPUT
% data:         n°Obs x n°Var
% subvar:       1 x n°Var (cell)
% color:        n°Var x 3


% OUTPUT
% h:        figure properties
% coord:    coord where the density distributions lay

nVar = size(data,2);

w = 0.03;          % scatterplot width
jit_base = 0;
base = 1;
coord(1) = base;

memo = [];
for n=1:nVar
    curr_data = data{n};
    [f(:,n), y(:,n)] = ksdensity(curr_data);
    means(n) = mean(curr_data);
    perc(n,:) = quantile(curr_data, [0.25 0.75]);
    if length(unique(curr_data))==1
        memo = [memo, n];
    end

end
space = max(max(f))/3;
space_jit = max(max(f))/6;

figure;
for n=1:nVar
    if ismember(n, memo)
        h{3,n} = scatter(base,means(n),50,'MarkerFaceColor',col(2,:),...
            'MarkerEdgeColor',col(2,:));
        
    else
        
        % half violin
        vec1L = cat(1, base-f(:,n), base*ones(length(f(:,n)),1));   % x coord
        vec2L = cat(1, y(:,n), y(:,n));                             % y coord
        h{1,n} = fill(vec1L(1:100), vec2L(1:100), color{1});
        h{1,n}.LineWidth = 1;
        h{1,n}.EdgeColor = color{1};
        hold on
        
        % mean and line plot
        h{3,n} = scatter(base,means(n),50,'MarkerFaceColor',color{2},...
            'MarkerEdgeColor',color{2});
        hold on
        h{4,n} = line([base base], [perc(n,1) means(n)],...
            'Color', color{2}, 'LineWidth', 1);
        hold on
        h{5,n} = line([base base], [means(n) perc(n,2)],...
            'Color', color{2}, 'LineWidth', 1);
        
        % scatterplot
        jit = (rand(size(data{n})))*w;    % jitter (y-coord for scatterplot)
        jit = jit + base + space_jit;
        curr_subvar = subvar(:,n);
        tj = 1;
        for k=1:length(curr_subvar)
            if not(isempty(curr_subvar{k}))
                h{2,n} = scatter(jit(tj:tj+length(curr_subvar{k})-1), curr_subvar{k});
                h{2,n}.SizeData = 3;
                h{2,n}.MarkerFaceColor = color{3}(k,:);
                h{2,n}.MarkerEdgeColor = 'none';
%                 h{2,n}.MarkerFaceAlpha = 0.6;
                tj = tj + length(curr_subvar{k});
            end
        end
        
    end
    
    if n<nVar
        base = base + w + space_jit + space + max(f(:,n+1));
        coord(n+1) = base;
    end
end
