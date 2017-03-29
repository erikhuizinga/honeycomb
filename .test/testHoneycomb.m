%% Test honeycomb

clear
close all


%% Test basic input arguments
% Set number of elements
n = 1000;

% Generate x
x = randn(n, 1);
x = [-5; -5; 5; 5; x];
x = x(1 : n);

% Generate y
y = randn(n, 1);
y = [-5; 5; -5; 5; y];
y = y(1 : n);


% Test honeycomb
figure
honeycomb(x, y)
colorbar


%TODO
% - Test data on edges: 6 edges between hexagons, outer edges of entire grid,
%   corners of hexagons