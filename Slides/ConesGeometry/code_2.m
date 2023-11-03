clc; close all; clear;

Count1 = 30;
Count2 = 30;

n = 3;
A = randn(3, 3);
A = A / norm(A);
iA = pinv(A);

c = randn(3, 1);

cc = iA'*c;

T = null(cc');

h = 1;
pn = h * cc / (norm(cc)*norm(cc));
gamma = norm(pn);

mu = sqrt(h^2 - gamma^2);

% xt = T*randn(2, 1);
% xt = mu*xt / norm(xt);

dphi = 2*pi / (Count2 - 1);

for j = 1:Count2
    
    phi = dphi * (j - 1);
    
    pt = mu*T*[cos(phi), sin(phi); -sin(phi), cos(phi)]*[1; 0];
    p = pt + pn;
    
    x = iA*p;
    
    [norm(p) - (cc'*p), ...
    norm(A*x) - (c'*x)]
    
%     res.X(i, j) = x(1);
%     res.Y(i, j) = x(2);
%     res.Z(i, j) = x(3);
end

% figure('Color', 'w');
% surf(res.X, res.Y, res.Z, 'FaceAlpha', 1.0, 'EdgeAlpha', 0.2); hold on;
% % surf(res.X, res.Y, res.Z, 'FaceColor', [0.2, 0.7, 0.3], 'FaceAlpha', 1.0, 'EdgeAlpha', 0.2);
% 
% plot3([0; 0.5], [0; 0], [0; 0], 'k', 'LineWidth', 0.5); text(0.6, 0, 0, 'x_1', 'FontSize', 13);
% plot3([0; 0], [0; 0.5], [0; 0], 'k', 'LineWidth', 0.5); text(0, 0.6, 0, 'x_2', 'FontSize', 13);
% plot3([0; 0], [0; 0], [0; 1.3], 'k', 'LineWidth', 0.5); text(0, 0, 1.4, 'J', 'FontSize', 13);
% 
% axis equal;
% axis off;

% xlabel('x_1');
% ylabel('x_2');
% zlabel('J');