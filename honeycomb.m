function varargout = honeycomb(varargin)
%HONEYCOMB Plot a bivariate histogram using hexagons a.k.a. honeycomb plot.
%
%   HONEYCOMB(X, Y) plots a bivariate histogram of X and Y using hexagonal
%   bins. X and Y can be arrays of any shape, but they must have the same
%   size. HONEYCOMB determines the bin edges using the same automatic
%   binning algorithm as histogram2, which returns uniform bins of an area
%   chosen to cover the range of values in X and Y and reveal the shape of
%   the underlying distribution.
%
%   HONEYCOMB(X, Y, NBINS), where NBINS is a scalar or 2-element vector,
%   specifies the number of bins to use. A scalar specifies the same number
%   of bins in each dimension, whereas a 2-element vector specifies a
%   different number of bins for the X and Y dimensions with the two
%   elements respectively.
%
%   P = HONEYCOMB(...) also returns a Patch object. Use this to inspect and
%   adjust properties of the hexagons.

% Copyright (c) 2017, Erik Huizinga


%% Validate input arguments
% Set defaults
defaults.numBins = [];
defaults.Debug = false;

% Keep track of input argument index
argIndex = 1;

% Instantiate an input parser
parser = inputParser;

% Add x and y arguments
parser.addRequired('x', @(x) validateData(x, mfilename, 'x', argIndex));
argIndex = argIndex + 1;
parser.addRequired('y', @(x) validateData(x, mfilename, 'y', argIndex));

% Add optional argument numBins: the number of hexagonal bins in both
% directions or, if specified as a two-element vector, the number of
% hexagonal bins in horizontal and vertical direction
argIndex = argIndex + 1;
parser.addOptional('numBins', defaults.numBins, ...
    @(x) validateNumBins(x, mfilename, argIndex));

% Add undocumented Name-Value pair argument to create additional debugging
% plots
parser.addParameter('Debug', defaults.Debug, @islogical)


% Parse input arguments
parser.parse(varargin{:});

% Get variables from input parser
struct2variables(parser.Results);
%#ok<*NODEF>
isDebug = Debug;

% Assert x and y have the same size
assert(all(size(x) == size(y)), ...
    'honeycomb:dimensionMismatch', ...
    'x and y must have the same dimensions')


%% Parse input arguments
% Parse data, exclude non-finite elements
isFinite = isfinite(x) & isfinite(y);
x(~isFinite) = [];
y(~isFinite) = [];

% Handle special number of elements
if isempty(x)
    warning( ...
        'honeycomb:noFiniteData', ...
        ['all data in x and / or y are non-finite; ' ...
        'plotting one hexagon only'])
end


% Parse number of bins
if isscalar(numBins)
    % Set both horizontal and vertical number of bins to numBins
    numXBins = numBins;
    numYBins = numBins;
    
    % Set a boolean to make it possible to draw one hexagon only
    isOneXBin = numBins == 1;
    isOneYBin = isOneXBin;
    
elseif isvector(numBins)
    % Get horizontal and vertical number of bins from numBins elements
    numXBins = numBins(1);
    numYBins = numBins(2);
    
    % Set a boolean to make it possible to draw one hexagon only
    isOneXBin = numXBins == 1;
    isOneYBin = numYBins == 1;
    
elseif numel(x) <= 1 && isempty(numBins)
    % Ensure the one hexagon is drawn
    numBins = 1;
    numXBins = 1;
    numYBins = 1;
    
    % Set a boolean to make it possible to draw one hexagon only
    isOneXBin = true;
    isOneYBin = true;
end
% At this point numXBins and numYBins are unknown if numBins is empty


%% Prepare number of hexagons
% Get square bin edges
if isempty(numBins)
    % Let the histcounts2 algorithm do the automatic binning
    [~, xEdges, yEdges] = histcounts2(x, y);
    
    % Determine the number of bins in either direction
    numXBins = numel(xEdges) - 1;
    numYBins = numel(yEdges) - mod(numel(yEdges), 2);
    
    % Set a boolean to make it possible to draw one hexagon only
    isOneXBin = numXBins == 1;
    isOneYBin = numYBins == 1;
    
else
    % Let the histcounts2 algorithm do the automatic binning
    % Note that the number of bins is decremented by one if greater than
    % one, because the number of hexagonal bins is one more than the number
    % of square bins
    [~, xEdges, yEdges] = ...
        histcounts2(x, y, [numXBins - ~isOneXBin, numYBins - ~isOneYBin]);
end
% At this point numXBins and numYBins are known


% Get range of bins, i.e., the rectangular area in which to position the
% hexagonal bins
xLim = xEdges([1, end]);
yLim = yEdges([1, end]);


% Calculate number of hexagon radii in data limits
numXRadii = 3 / 2 * (numXBins - 1) + 1;
if isOneXBin && isOneYBin
    numYRadii = sqrt(3);
elseif isOneXBin
    numYRadii = sqrt(3) * numYBins;
else
    numYRadii = sqrt(3) / 2 * (2 * numYBins - 1);
end

% Calculate radii
rx = diff(xLim) / numXRadii;
ry = diff(yLim) / numYRadii;

% Calculate distance between hexagon centers
dx = 3 / 2 * rx;
dy = sqrt(3) * ry;


%% Calculate hexagon centers
% Note: the area spanning the horizontal data range is rx / 2 left and
% right of the centers
if isOneXBin
    xCenters = mean(xLim);
else
    xCenters = (xLim(1) + rx/2) : dx : (xLim(2) - rx/2);
end

% Set the vertical shift distance for every second hexagon column
if ~isOneXBin
    sy = dy / 2;
end

% Note: the area spanning the vertical data range is dy / 2 above and below
% the centers
if isOneXBin && ~isOneYBin
    % Space the centers such that the top and bottom hexagon vertices span
    % the entire range
    yCenters = linspace(yLim(1) + dy/2, yLim(2) - dy/2, numYBins);
    
elseif isOneYBin
    if isOneXBin
        % Set the center to the middle of the data range
        yCenters = mean(yLim);
        
    else
        % Set the center to the middle, but allow for the shift
        yCenters = mean(yLim) - sy/2;
    end
else
    % Set the default vertical centers
    yCenters = linspace(yLim(1), yLim(2) - sy, numYBins);
end


% Initialize a grid of hexagon centers
[XCenters, YCenters] = meshgrid(xCenters, yCenters);

% Shift y centers on evey second horizontal bin
if ~isOneXBin
    isShifted = 1 : 2 : numXBins;
    YCenters(:, isShifted) = YCenters(:, isShifted) + dy/2;
end


% Vectorize centers
XCenters = XCenters(:);
YCenters = YCenters(:);

% Sort center coordinates from top to bottom and right to left
Centers = flipud(sortrows([YCenters, XCenters]));
YCenters = Centers(:, 1);
XCenters = Centers(:, 2);


%% Calculate hexagon vertices
% Set vertex angles from center
theta = (1/6 : 1/6 : 1) * 360;

% Calculate vertex coordinates
XVertices = bsxfun(@plus, XCenters, rx * cosd(theta)).';
YVertices = bsxfun(@plus, YCenters, ry * sind(theta)).';


%% Determine bin counts
% Get number of hexagons
numHexagons = numel(XCenters);

% Preallocate bin counts
counts = zeros(numHexagons, 1);

% Preallocate loop counter
n = 1;

% Copy data
xx = x;
yy = y;

% Loop over hexagons
while n <= numHexagons && ~isempty(xx)
    % Determine bin count
    isInBin = inpolygon(xx, yy, XVertices(:, n), YVertices(:, n));
    counts(n) = nnz(isInBin);
    
    % Clear bin data to prevent double counts on edges
    xx(isInBin) = [];
    yy(isInBin) = [];
    
    % Increment counter
    n = n + 1;
end

% Determine where no data was counted
if isOneXBin && isOneYBin
    isIncluded = true;
else
    isIncluded = counts > 0;
end


%% Draw hexagons
% Prepare axes for a new plot
ax = newplot;

% Draw hexagons
p = patch( ...
    XVertices(:, isIncluded), ...
    YVertices(:, isIncluded), ...
    counts(isIncluded));


%% Customize axes
if ~ishold(ax)
    % Set axes properties
    axis(ax, 'tight')
    ax.Box = 'on';
end


%% Create debug plots
if isDebug
    hold on
    
    % Scatter data
    plot(x(:), y(:), 'or')
    
    % Draw data limits
    plot( ...
        [xLim(1), xLim, flip(xLim)], [yLim, flip(yLim), yLim(1)], ...
        'rs:', 'MarkerFaceColor', 'r')
    
    % Scatter hexagon centers and vertices
    plot(XCenters, YCenters, 'ks')
    plot(XVertices(:), YVertices(:), 'kv')
end


%% Set output
if nargout
    varargout = {p};
end
end