%% Test honeycomb

clear
close all


%% Test discretize
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

% Set number of hexagons in x and y direction
numBins = 10;

% Plot bivariate histogram viewed from the top
figure
hold on
histogram2(x, y, numBins, 'DisplayStyle', 'tile')
plot(x, y, 'ro')

% Test honeycomb
figure
honeycomb(x, y, numBins)
colorbar

%TODO
% - Test data on edges: 6 edges between hexagons, outer edges of entire grid,
%   corners of hexagons