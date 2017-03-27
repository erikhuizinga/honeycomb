function honeycomb(x, y, numBins)
% HONEYCOMB


%% Discretize data into bins
% Get number of bins in both directions
numXBins = numBins;
numYBins = numBins;

% Get square bin edges
% The number of bins is decremented by one for the discretize function,
% because the number of hexagonal bins is one more than the number of
% square bins
[~, xEdges] = discretize(x, numXBins - 1);
[~, yEdges] = discretize(y, numYBins - 1);

% Get range of bins
xLim = xEdges([1, end]);
yLim = yEdges([1, end]);

% Calculate number of hexagon radii in range
numRx = 3 / 2 * (numXBins - 1) + 1;
numRy = 2 * numYBins - 1;

% Set hexagon radius
rx = diff(xLim) / numRx;
ry = diff(yLim) / numRy;

% Set the distance between hexagon centers
dx = 3 / 2 * rx;
dy = 3 / 2 * ry;


% Calculate hexagon centers
xCenter = linspace(xLim(1) + rx/2, xLim(2) - rx/2, numXBins);
yCenter = yEdges;


%% Initialize a grid of hexagon centers
[XCenter, YCenter] = meshgrid(xCenter, yCenter);

% Shift y centers
isShifted = 1 : 2 : numYBins;
YCenter(:, isShifted) = YCenter(:, isShifted) + dy/2;

% Vectorize centers
XCenter = XCenter(:);
YCenter = YCenter(:);


%TODO delete me
figure
plot(XCenter, YCenter, 'r^')


%% Initialize hexagon vertices
theta = linspace(pi / 3, 2 * pi, 6);
XVertices = bsxfun(@plus, XCenter, rx * cos(theta));
YVertices = bsxfun(@plus, YCenter, ry * sin(theta));


%TODO delete me
hold on
plot(XVertices(:), YVertices(:), 'bo')

end