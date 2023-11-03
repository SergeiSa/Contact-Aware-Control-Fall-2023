clc; close all; clear;

Count1 = 30;
Count2 = 30;

res.X = zeros(Count1, Count2);
res.Y = zeros(Count1, Count2);
res.Z = zeros(Count1, Count2);

dr = 1 / Count1;
dphi = 2*pi / (Count2 - 1);

for i = 1:Count1
for j = 1:Count2
    r = dr * (i - 1);
    phi = dphi * (j - 1);
    
    x1 = r*cos(phi);
    x2 = -r*sin(phi);
    z = norm([x1, x2]);
    
    res.X(i, j) = x1;
    res.Y(i, j) = x2;
    res.Z(i, j) = z;
end
end

figure('Color', 'w');
surf(res.X, res.Y, res.Z, 'FaceAlpha', 1.0, 'EdgeAlpha', 0.2); hold on;
% surf(res.X, res.Y, res.Z, 'FaceColor', [0.2, 0.7, 0.3], 'FaceAlpha', 1.0, 'EdgeAlpha', 0.2);

plot3([0; 0.5], [0; 0], [0; 0], 'k', 'LineWidth', 0.5); text(0.6, 0, 0, 'x_1', 'FontSize', 13);
plot3([0; 0], [0; 0.5], [0; 0], 'k', 'LineWidth', 0.5); text(0, 0.6, 0, 'x_2', 'FontSize', 13);
plot3([0; 0], [0; 0], [0; 1.3], 'k', 'LineWidth', 0.5); text(0, 0, 1.4, 'J', 'FontSize', 13);

axis equal;
axis off;

% xlabel('x_1');
% ylabel('x_2');
% zlabel('J');