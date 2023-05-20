clc
clear

% Define the range of x and y values
x = linspace(-5, 5, 100);
y = linspace(-5, 5, 100);

% Create a grid of points
[X, Y] = meshgrid(x, y);

% Define parameters for the paraboloid equation
a = 1;  % Coefficient of x^2
b = 1;  % Coefficient of y^2
c = 1;  % Coefficient of z

% Calculate the z-values based on the paraboloid equation
Z = a*X.^2 + b*Y.^2 + c;

% Create the surface plot
figure;
surf(X, Y, Z);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Paraboloid Surface');
