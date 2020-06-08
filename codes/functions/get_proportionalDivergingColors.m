function [propColors, limCmap] = get_proportionalDivergingColors(mat, colors)

% given a colormap and a vector in which in every row there's a value of a
% certain variable, this function return a set of colors proportional to
% such values.

% INPUT:    mat = nObs*nVar
%           colors = nColColormap*3
%
% OUTPUT:   propColors = nObs*3

if nargin<2
    colors = cbrewer('div','PiYG',201,'pchip');
end

curr_max = zeros(size(mat,2),1);
curr_min = zeros(size(mat,2),1);
for i=1:size(mat,2)
    curr_max(i) = max(mat(:,i));
    curr_min(i) = min(mat(:,i));
end
clear n

diffPos = max(curr_max);
diffNeg = min(curr_min);

limCmap = max([diffPos abs(diffNeg)]);

ind_zero = round(size(colors,1)/2);

propColors = zeros(size(mat,1), 3, size(mat,2));
for i=1:size(mat,2)
    for j=1:size(mat,1)
        if mat(j,i)>0
            tmp = mat(j,i)*(floor(size(colors,1))/2)/limCmap;
            if mat(j,i) == limCmap
                ind = size(colors,1);
            else
                ind = round(ind_zero+tmp);
            end
            propColors(j,:,i) = colors(ind,:);
        elseif mat(j,i)<0
            tmp = abs(mat(j,i))*(floor(size(colors,1))/2)/limCmap;
            if abs(mat(j,i)) == limCmap
                ind = 1;
            else
                ind = round(ind_zero-tmp);
                if ind==0
                    ind = 1;
                end
            end
            propColors(j,:,i) = colors(ind,:);
        elseif mat(j,i)==0
            propColors(j,:,i) = colors(ind_zero,:);
        end        
    end
    clear j
end
clear i