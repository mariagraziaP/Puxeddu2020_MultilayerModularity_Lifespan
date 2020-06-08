function [propColors, limCmap] = get_proportionalColors(mat, colors)

% given a colormap and a matrix in which in every row there's a value of a
% certain variable, this function return a set of colors proportional to
% such values.

% INPUT:    mat = nObs*nVar
%           colors = nColColormap*3
%
% OUTPUT:   propColors = nObs*3

if nargin<2
    colors = cbrewer('seq','Reds',1000,'pchip');
end

curr_max = zeros(size(mat,2),1);
curr_min = zeros(size(mat,2),1);
for i=1:size(mat,2)
    curr_max(i) = max(mat(:,i));
    curr_min(i) = min(mat(:,i));
end
clear n
maxV = max(curr_max);
minV = min(curr_min);
diffV = maxV-minV;

limCmap = [minV maxV];

propColors = zeros(size(mat,1), 3, size(mat,2));
for i=1:size(mat,2)
    for j=1:size(mat,1)
        ind = round(((mat(j,i)-minV)*size(colors,1))/diffV); 
        if ind==0
            ind = 1;
        end
        propColors(j,:,i) = colors(ind,:);
    end
    clear j
end
clear i