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

% Calculate number of hexagon radii in data limits and calculate distance
% between hexagon centers
numXRadii = 3 / 2 * (numXBins - 1) + 1;
rx = diff(xLim) / numXRadii;
dx = 3 / 2 * rx;

numYRadii = sqrt(3) / 2 * (2 * numYBins - 1);
ry = diff(yLim) / numYRadii;
dy = sqrt(3) * ry;


% Calculate hexagon centers, such that hexagons cover entire data range
% from x,y limit 1 to x,y limit 2
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


%TODO delete me
figure
plot(XCenter, YCenter, 'r^')


%% Initialize hexagon vertices
theta = (1/6 : 1/6 : 1) * 360;
XVertices = bsxfun(@plus, XCenter, rx * cosd(theta));
YVertices = bsxfun(@plus, YCenter, ry * sind(theta));


%TODO delete me
hold on
plot(XVertices(:), YVertices(:), 'bo')


end