clc; close all; clear;

n = randn(3, 1);
n = n / norm(n);
% n = [0;0;1];
mu = 0.6;

T = null(n');
P = eye(3) - n*n';

Count = 10;
dphi = 2*pi / (Count - 1);    

Res.V = zeros(3, Count);
Res.normals = zeros(3, Count-1);

Res.X = zeros(2, Count);
Res.Y = zeros(2, Count);
Res.Z = zeros(2, Count);

%generation

for j = 1:Count
        
        phi = dphi * (j - 1);
        point = n + mu * T*[cos(phi); -sin(phi)];
        
        Res.V(:, j) = point;
        Res.X(1, j) = point(1); Res.X(2, j) = 0;
        Res.Y(1, j) = point(2); Res.Y(2, j) = 0;
        Res.Z(1, j) = point(3); Res.Z(2, j) = 0;
end


for j = 1:(Count-1)
    
    index1 = j;
    index2 = j+1;
    
    point1 = Res.V(:, index1);
    point2 = Res.V(:, index2);
    
    a = cross(point1, point2);
    a = a / norm(a);
    
%     a = a * sign(dot(P*a, P*(point1+point2)/2));
    a = a * sign(a' * P * (point1+point2));
    
    Res.normals(:, j) = a;
end
A = Res.normals';
b = zeros(Count-1, 1);

%test

for i = 1:10000
    
    x = randn(3, 1);
    
    check_SOC = (norm(P*x) <= mu*dot(n, x));
    if A*x <= b
        check_H = true;
    else
        check_H = false;
    end
    
    if check_H &&(~check_SOC)
        error('inner approximation failed');
    end
    if check_SOC &&(~check_H)
        disp('approximation discrepancy detected');
    end
end

%graphics

figure('Color', 'w');
for j = 1:Count
    temp_point = Res.V(:, j);
    
    plot3([0, temp_point(1)], [0, temp_point(2)], [0, temp_point(3)], 'LineWidth', 0.5, 'Color', 'k'); hold on;
end
for j = 1:(Count-1)
    index1 = j;
    index2 = j+1;
    
    temp_point1 = Res.V(:, index1);
    temp_point2 = Res.V(:, index2);
    temp_point = 0.5 * (temp_point1 + temp_point2);
    temp_norm = 0.1*Res.normals(:, j) + temp_point;
    
    plot3([temp_point(1), temp_norm(1)], ...
          [temp_point(2), temp_norm(2)], ...
          [temp_point(3), temp_norm(3)], 'LineWidth', 1, 'Color', 'r'); hold on;
end
plot3(Res.V(1, :), Res.V(2, :), Res.V(3, :), 'o');
surf(Res.X, Res.Y, Res.Z, 'FaceAlpha', 0.5);

grid on;
axis equal;