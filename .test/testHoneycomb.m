%% Test honeycomb
%TODO
% - Test data on edges: 6 edges between hexagons, outer edges of entire grid,
%   corners of hexagons


%% Initialize
clear
close all


%% Test basic input arguments
% Set number of elements
n = 50;  % Will be squared
c = 3;

% Generate x
x = randn(n);
preX = c * [-1; -1; 1; 1];
x(1 : numel(preX)) = preX;

% Generate y
y = randn(size(x));
preY = c * [-1; 1; -1; 1];
y(1 : numel(preY)) = preY;


% Test honeycomb
figure
honeycomb(x, y)
title 'Default'
colorbar


%% Test data validation and parsing
% Generate unplottable data
xNaN = NaN(size(x));
xNaN(end) = 1;
yNaN = y;
yNaN(end) = NaN;

% Test unplottable data
figure
honeycomb(xNaN, yNaN, 'Debug', true)
title 'Only NaN'
colorbar


% Test empty data
% This must return a valid patch object, e.g., compare to
% h = scatter([], [])
figure
p = honeycomb([], [], 'Debug', true);
title 'Empty data'
colorbar
assert(isa(p, 'matlab.graphics.primitive.Patch'), ...
    'output argument of honeycomb must be a Patch object')


% Test singleton data
figure
honeycomb(-1, pi, 'Debug', true)
title 'Singleton data'
colorbar


% Test singleton data with specifed number of bins
figure
honeycomb(-pi, 1, 1, 'Debug', true)
title 'Singleton data in one bin'
colorbar

figure
honeycomb(exp(1), sqrt(2), 2, 'Debug', true)
title 'Singleton data in 2 bins'
colorbar

figure
honeycomb(exp(-2), -sqrt(3), [3, 2], 'Debug', true)
title 'Singleton data in [3, 2] bins'
colorbar


% Test invalid data
isExceptionThrown = false;
try
    figure
    honeycomb(1i, 2, 'Debug', true)
catch
    title 'Invalid data, empty plot'
    isExceptionThrown = true;
end
assert(isExceptionThrown, 'an exception must be thrown')


%% Test validation and parsing of number of bins
% Test drawing one hexagonal bin
figure
honeycomb(x, y, 1, 'Debug', true)
title 'One bin'
colorbar

% Test drawing lots of bins
figure
honeycomb(x, y, 100, 'Debug', true)
title '100 bins'
colorbar

% Test drawing one horizontal hexagonal bin against many vertical bins
figure
honeycomb(x, y, [1, 10], 'Debug', true)
title '[1, 10] bins'
colorbar

% Test drawing one vertical hexagonal bin against many horizontal bins
figure
honeycomb(x, y, [10, 1], 'Debug', true)
title '[10, 1] bins'
colorbar
