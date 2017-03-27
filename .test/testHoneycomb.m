%% Test honeycomb

clear
close all


%% Test discretize
% Set number of elements
n = 100;

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

% Test honeycomb
honeycomb(x, y, numBins)