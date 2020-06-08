function [l, b] = plot_BoundedLines(matrix,colors)

% This function returns a graph plot of the means of the variables with
% the respective confidence intervals, shown through a shaded area 

% INPUT:    matrix = [n°xpoints, n°obs, n°var]
%           colors = [n°var, 3, 2]
%
% OUTPUT:   l = properties of lines
%           b = properties of shaded area

nVar = size(matrix,3);
nPoints = size(matrix,1);

if nargin<2
    col = cbrewer('qual','Paired',nVar*2,'pchip');
    i = 1;
    for n=1:nVar
        colors(n,:,1) = col(i+1,:);
        colors(n,:,2) = col(i,:);
        i = i+2;
    end
    clear n
end

[mean_var, bound_var] = findCIforBoundedLine(matrix);

l = cell(nVar,1);
b = cell(nVar,1);

x = 1:nPoints;
xbound = [x, fliplr(x)];

% figure;
for n=1:nVar
    y1 = mean_var(:,n)-bound_var(:,1,n);
    y1 = y1';
    y2 = mean_var(:,n)+bound_var(:,2,n);
    y2 = y2';
    ybound = [y1, fliplr(y2)];
    b{n} = fill(xbound, ybound, colors(n,:,2));
    b{n}.EdgeColor = 'none';
    b{n}.FaceAlpha = 0.5;
    hold on
    l{n} = plot(1:nPoints,mean_var(:,n),'Color',colors(n,:,1));
    hold on
    clear y1 y2 ybound
end
xlim([1 nPoints])