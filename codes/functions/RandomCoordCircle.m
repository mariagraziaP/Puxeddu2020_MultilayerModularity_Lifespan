function [x,y] = RandomCoordCircle(cx,cy,radius,points)

% function that find random coordinates within a circle

% INPUT: 
% cx,cy = centre coordinates
% radius = circle radius
% points = number of random coordinates desired

% OUTPUT:
% x,y = a number equal to 'points' of random coordinates within the circle of
%       radius equal to 'radius' and centre in [cx,cy]

angle1 = 0;
angle2 = 2*pi;
t = (angle2 - angle1) * rand(points,1) + angle1;
r = radius*sqrt(rand(points,1));
x = cx + r.*cos(t);
y = cy + r.*sin(t);
end




% % Create a random set of coordinates in a circle.
% % First define parameters that define the number of points and the circle.
% n = 5000;
% R = 20;
% x0 = 50; % Center of the circle in the x direction.
% y0 = 90; % Center of the circle in the y direction.
% 
% % Now create the set of points.
% % For a full circle, use 0 and 2*pi.
% angle1 = 0;
% angle2 = 2*pi;
% % For a sector, use partial angles.
% % angle1 = pi/4;
% % angle2 = 3*pi/4;
% t = (angle2 - angle1) * rand(n,1) + angle1;
% r = R*sqrt(rand(n,1));
% x = x0 + r.*cos(t);
% y = y0 + r.*sin(t);
% r = R*sqrt(rand(n,1));
% x = x0 + r.*cos(t);
% y = y0 + r.*sin(t);
% 
% % Now display our random set of points in a figure.
% plot(x,y, '.', 'MarkerSize', 5)
% axis square;
% grid on;
% % Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% fontSize = 30;
% xlabel('X', 'FontSize', fontSize);
% ylabel('Y', 'FontSize', fontSize);
% title('Random Locations Within a Circle', 'FontSize', fontSize);