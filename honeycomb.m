function varargout = honeycomb(varargin)
% HONEYCOMB


%% Validate input arguments
% Set defaults
defaults.numBins = [];

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

% Parse input arguments
parser.parse(varargin{:});

% Get variables from input parser
struct2variables(parser.Results);
%#ok<*NODEF>


%% Parse input arguments
% Parse number of bins
if isscalar(numBins)
    numXBins = numBins;
    numYBins = numBins;
elseif isvector(numBins)
    numXBins = numBins(1);
    numYBins = numBins(2);
end


%% Discretize data into bins
% Get square bin edges
% The number of bins is decremented by one for the discretize function,
% because the number of hexagonal bins is one more than the number of
% square bins
if isempty(numBins)
    [~, xEdges, yEdges] = histcounts2(x, y);
    numXBins = numel(xEdges) - 1;
    numYBins = numel(yEdges) - mod(numel(yEdges), 2);
else
    [~, xEdges, yEdges] = histcounts2(x, y, [numXBins - 1, numYBins - 1]);
end

% Get range of bins, i.e., the rectangular area in which to position the
% hexagonal bins 
xLim = xEdges([1, end]);
yLim = yEdges([1, end]);


% Calculate number of hexagon radii in data limits and calculate distance
% between hexagon centers
numXRadii = 3 / 2 * (numXBins - 1) + 1;
numYRadii = sqrt(3) / 2 * (2 * numYBins - 1);

% Calculate radii
rx = diff(xLim) / numXRadii;
ry = diff(yLim) / numYRadii;

% Calculate distance between hexagon centers
dx = 3 / 2 * rx;
dy = sqrt(3) * ry;


% Calculate hexagon centers, such that hexagons cover entire data range
xCenter = (xLim(1) + rx/2) : dx : (xLim(2) - rx/2);
yCenter = yLim(1) : dy : (yLim(2) - dy/2);


%% Initialize a grid of hexagon centers
[XCenter, YCenter] = meshgrid(xCenter, yCenter);

% Shift y centers
isShifted = 1 : 2 : numYBins;
YCenter(:, isShifted) = YCenter(:, isShifted) + dy/2;

% Vectorize centers
XCenter = XCenter(:);
YCenter = YCenter(:);

% Sort center coordinates from top to bottom and right to left
Centers = flipud(sortrows([YCenter, XCenter]));
YCenter = Centers(:, 1);
XCenter = Centers(:, 2);


%% Initialize hexagon vertices
theta = (1/6 : 1/6 : 1) * 360;
XVertices = bsxfun(@plus, XCenter, rx * cosd(theta)).';
YVertices = bsxfun(@plus, YCenter, ry * sind(theta)).';


%% Determine bin counts
% Get number of hexagons
numHexagons = numel(XCenter);

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


%% Draw hexagons
% Prepare axes for a new plot
ax = newplot;

% Determine where no data was counted
isIncluded = counts > 0;

% Draw hexagons
p = patch( ...
    XVertices(:, isIncluded), ...
    YVertices(:, isIncluded), ...
    counts(isIncluded));


% Set axes limits
if ax.XTickMode == 'auto'
    xlim(xLim)
end
if ax.YTickMode == 'auto'
    ylim(yLim)
end

%% %TODO debug plots
hold on

% Scatter data
plot(x(:), y(:), 'or')

% Draw data limits
plot([xLim(1), xLim, flip(xLim)], [yLim, flip(yLim), yLim(1)], 'r:')

% Scatter hexagon centers and vertices
plot(XCenter, YCenter, 'ks')
plot(XVertices(:), YVertices(:), 'kv')


%% Set output
if nargout
    varargout = {p};
end
end