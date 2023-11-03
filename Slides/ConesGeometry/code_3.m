clc; close all; clear;

Count1 = 30;
Count2 = 30;

n = 3;

A = randn(3, 3);
A = A / norm(A);
iA = pinv(A);

c = randn(3, 1);
c = c / norm(c);
cc = iA'*c;

T = null(cc');

res.X = zeros(Count1, Count2);
res.Y = zeros(Count1, Count2);
res.Z = zeros(Count1, Count2);

res.line = zeros(Count1, 3);

dh = 1 / Count1;
dphi = 2*pi / (Count2 - 1);

for i = 1:Count1
    
    h = dh * (i - 1);
    
    pn = h * cc / (norm(cc)*norm(cc));
    gamma = norm(pn);
    
    mu = sqrt(h^2 - gamma^2);
    
    for j = 1:Count2
        
        phi = dphi * (j - 1);
        
        pt = mu * T*[cos(phi), sin(phi); -sin(phi), cos(phi)] * [1; 0];
        p = pt + pn;
        
        x = iA*p;
        x0 = iA*pn;
        
        if (abs(norm(p) - (cc'*p)) > 0.001) || ...
            abs(norm(A*x) - (c'*x)) > 0.001
            warning('discrepancy!')
        end
        
        res.X(i, j) = x(1);
        res.Y(i, j) = x(2);
        res.Z(i, j) = x(3);
        
        res.line(i, :) = x0;
    end
end

Points = [res.X(:), res.Y(:), res.Z(:)]';
APoints = A*Points;
normAPoints = vecnorm(APoints);
maxnormAPoints = max(normAPoints);

res.C = zeros(Count1, Count2, 3);
for i = 1:Count1
for j = 1:Count2
    x = [res.X(i, j); res.Y(i, j); res.Z(i, j)];
    normAx = norm(A*x);
    res.C(i, j, 1) = 0.1;
    res.C(i, j, 2) = normAx / maxnormAPoints;
    res.C(i, j, 3) = 0.8*normAx / maxnormAPoints;
end
end

figure('Color', 'w');
surf(res.X, res.Y, res.Z, res.C, 'FaceAlpha', 1.0, 'EdgeAlpha', 0.2); hold on;
% surf(res.X, res.Y, res.Z, 'FaceColor', [0.2, 0.7, 0.3], 'FaceAlpha', 1.0, 'EdgeAlpha', 0.2);

plot3([0; 0.5], [0; 0], [0; 0], 'k', 'LineWidth', 0.5); text(0.6, 0, 0, 'x_1', 'FontSize', 13);
plot3([0; 0], [0; 0.5], [0; 0], 'k', 'LineWidth', 0.5); text(0, 0.6, 0, 'x_2', 'FontSize', 13);
plot3([0; 0], [0; 0], [0; 1.3], 'k', 'LineWidth', 0.5); text(0, 0, 1.4, 'x_3', 'FontSize', 13);

plot3([0; c(1)], [0; c(2)], [0; c(3)], 'r', 'LineWidth', 1.0); text(0, 0, 1.4, 'c', 'FontSize', 13);

plot3(res.line(:, 1), res.line(:, 2), res.line(:, 3), '--', 'Color', [0.5, 0.1, 0.1], 'LineWidth', 1.0)

axis equal;
axis off;

% xlabel('x_1');
% ylabel('x_2');
% zlabel('J');