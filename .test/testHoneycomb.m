%% Test honeycomb

clear
close all


%% Test discretize
% Set number of elements
n = 30;

% Generate x
x = rand(n, 1);
x = [0; 0; 1; 1; x];
x = x(1 : n);

% Generate y
y = rand(n, 1);
y = [0; 1; 0; 1; y];
y = y(1 : n);

% Set number of hexagons in x and y direction
numBins = 3;

% Plot bivariate histogram viewed from the top
figure
hold on
histogram2(x, y, numBins, 'DisplayStyle', 'tile')
plot(x, y, 'ro', 'MarkerFaceColor', 'r')

% Test honeycomb
honeycomb(x, y, numBins)


%TODO
% - Test data on edges: 6 edges between hexagons, outer edges of entire grid,
%   corners of hexagons