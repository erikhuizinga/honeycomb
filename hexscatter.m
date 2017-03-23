function h = hexscatter(xData, yData, varargin)
% HEXSCATTER Create a scatter plot using hexagonal binning
%
%   HEXSCATTER(x, y) creates a scatter plot using hexagonal binning of the
%   data in x and y. The colour intensity of the hexagons indicates the
%   density of points in the hexagon. x and y must contain the same number
%   of elements and can be vectors or matrices, in which every element in x
%   belongs to the same element in y (using linear indexing in the case of
%   matrices). Any NaN or Inf values in x or y are ignored (with the
%   corresponding element in the other array). Only the real part of any
%   complex numbers in x or y will be used an a warning will be issued.
%
%   HEXSCATTER(..., Name, Value) sets options for the plot using Name-Value
%   pair arguments. The possible arguments are listed below.
%
%   h = HEXSCATTER(...) returns the handles to the created patch object,
%   which is used to implement this function.
%
%   Name-Value Pair Arguments
%   Specify optional comma-separated pairs of Name,Value arguments to
%   access various options. Name is the argument name and Value is the
%   corresponding value. Name must appear inside single quotes (' '). You
%   can specify several name and value pair arguments in any order as
%   Name1, Value1, ..., NameN, ValueN.
%   Example: 'XLim', [0, 1], 'YLim', [1, 10]
%
%   'XLim': x-axis limits
%   [min(x(:)), max(x(:))] (default) | two-element numeric vector
%   The x-axis limits specified as a two-element numeric vector. The first
%   and second element correspond to the left and right x-axis limit
%   respectively. The second element must be greater than the first.
%
%   'YLim': y-axis limits
%   [min(y(:)), max(y(:))] (default) | two-element numeric vector
%   The y-axis limits specified as a two-element numeric vector. The first
%   and second element correspond to the lower and upper y-axis limit
%   respectively. The second element must be greater than the first.
%
%   'Resolution': resolution of hexagons
%   [100, 100] (default) | numeric scalar | two-element numeric vector
%   The hexagonal grid resolution specified either as a scalar, or a
%   two-element vector. The resolution is the number of hexagons in
%   horizontal and vertical direction. If the specified value is a
%   two-element vector, the first element corresponds to the vertical
%   resolution and the second to the horizontal resolution. If the
%   specified value is scalar, then its value is used for both directions.
%
%   'DrawEdges': draw edges around hexagons
%   false (default) | true
%   Draw edges around the hexagons if the specified value is true.
%
%   'ShowZeros': show zero count hexagons
%   false (default) | true

% This function is based on hexscatter.m by Gordon Bean, source:
% https://raw.githubusercontent.com/brazilbean/bean-matlab-toolkit/
% ab1885c10f6bf96307fcbf73b6eea7657bf38da0/hexscatter.m
%
% Also available in the Bean Matlab Toolkit:
% https://github.com/brazilbean/bean-matlab-toolkit
%
% License of original:
% https://raw.githubusercontent.com/brazilbean/bean-matlab-toolkit/
% 1b761acf84330dca577ac0e60e938dc785989acc/LICENSE

% Copyright (c) 2017, Salman Mashayekh, Gordon Bean, Erik Huizinga
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the
%       distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


%% Parse input arguments
% Convert to vectors
xData = xData(:);
yData = yData(:);

% Keep only valid indices
lValid = isnan(xData) | isnan(yData) | isinf(xData) | isinf(yData);
xData = xData(~lValid);
yData = yData(~lValid);


% Parse variable input arguments
parameters = getDefaultParameters( ...
    varargin, ...
    'XLim', [min(xData(:)), max(xData(:))], ...
    'YLim', [min(yData(:)), max(yData(:))], ...
    'Resolution', [100, 100], ...
    'DrawEdges', false, ...
    'ShowZeros', false ...
    );


% Format resolution parameters
if length(parameters.Resolution) == 1
    parameters.Resolution = parameters.Resolution * [1, 1];
end
parameters.x.resolution = parameters.Resolution(2);
parameters.y.resolution = parameters.Resolution(1);


% Format limit parameters
parameters.x.limit = parameters.XLim;
parameters.y.limit = parameters.YLim;


%% Determine hexagon grid
% Determine hexagon centre coordinates
xCenter = linspace( ...
    parameters.x.limit(1), parameters.x.limit(2), parameters.x.resolution);
yCenter = linspace( ...
    parameters.y.limit(1), parameters.y.limit(2), parameters.y.resolution);

% Make X and Y grids of hexagon centres
[XCenter, YCenter] = meshgrid(xCenter, yCenter);

% Adjust every second hexagon centre in Y
dy = diff(yCenter([1, 2])) / 2;
iAdjust = 1 : 2 : numel(xCenter);
YCenter(:, iAdjust) = YCenter(:, iAdjust) + dy;
yCenter(iAdjust) = yCenter(iAdjust) + dy; %FIXME


%% Map points to boxes
% Determine horizontal and vertical box numbers of data
dx = diff(xCenter([1, 2]));
boxX = floor((xData - parameters.x.limit(1)) / dx) + 1;
boxY = floor((yData - parameters.y.limit(1)) / (2 * dy)) + 1;


% Limit box numbers to hexagonal grid range
%TODO Ignore data points outside grid: user requests grid of a limited
%     size, i.e., points outside it should not be counted in the grids
boxX = min(boxX, parameters.x.resolution);
boxX = max(boxX, 1);
boxY = min(boxY, parameters.x.resolution);
boxY = max(boxY, 1);


% Determine orientation
orientation = mod(boxX, 2) == 1;

% Shift data to relative position in boxs
boxXData = xData - xCenter(boxX).';
boxYData = yData - yCenter(boxY).';

% Which layer
layer = boxYData > dy;

% Convert to block B format
l2Flip = layer == orientation;
boxXData(l2Flip) = dx - boxXData(l2Flip);

boxYData(layer == 1) = boxYData(layer == 1) - dy;

%TODO
newplot
hold on
axis equal
plot(boxXData, boxYData, 'rx')

% Find closest corner
topright = ...
    boxXData.^2 + boxYData.^2 > (boxXData - dx).^2 + (boxYData - dy).^2;

%TODO
table(xData, yData, boxXData, boxYData, boxX, boxY, orientation, layer, topright)

%% Map corners back to bins
% Determine horizontal bin number
binX = boxX + ~(orientation == (layer == topright));
binX(binX > length(xCenter)) = length(xCenter);

% Determine vertical bin number
binY = boxY + (layer & topright);
binY(binY > length(yCenter)) = length(yCenter);

binX = sub2ind(size(XCenter), binY, binX);
% Replaced foox to conserve memory


%% Determine counts
counts = histc(binX, 1 : numel(XCenter));

%newplot; %TODO
xscale = dx * 2 / 3;
yscale = dy * 2 / sqrt(3);
theta = (0 : 1/3 : 2) * pi;
x = bsxfun(@plus, XCenter(:), cos(theta) * xscale).';
y = bsxfun(@plus, YCenter(:), sin(theta) * yscale).';

%TODO
plot(XCenter(:), YCenter(:), 'k^')
plot(x(:), y(:), 'ko')

% Determine patch edgeColor value
if parameters.DrawEdges
    patchEdgeColor = 'flat';
else
    patchEdgeColor = 'none';
end

plot(xData, yData, 'ro')


% Draw patches
if parameters.ShowZeros
    h = patch(x, y, counts, 'edgeColor', patchEdgeColor);
else
    jj = counts > 0;
    h = patch(x(:, jj), y(:, jj), counts(jj), ...
        'EdgeColor', patchEdgeColor, ...
        'FaceAlpha', .2);
end

xlim(parameters.XLim);
ylim(parameters.YLim);

if nargout == 0
    clear h;
end

%% Function: getDefaultParameters
% Gordon Bean, March 2012
% Copied from https://github.com/brazilbean/bean-matlab-toolkit
    function parameters = getDefaultParameters(parameters, varargin)
        if (iscell(parameters))
            parameters = getParameters(parameters{:});
        end
        defaults = getParameters(varargin{:});
        
        for f = fieldnames(defaults).'
            field = f{:};
            if (~isfield(parameters, field))
                parameters.(field) = defaults.(field);
            end
        end
    end

%% Function: get_params - return a struct of key-value pairs
% Gordon Bean, January 2012
%
% Usage
% params = get_params( ... )
%
% Used to parse key-value pairs in varargin - returns a struct.
% Converts all keys to lower case.
%
% Copied from https://github.com/brazilbean/bean-matlab-toolkit
    function parameters = getParameters(varargin)
        parameters = struct;
        
        nn = length(varargin);
        if (mod(nn, 2) ~= 0)
            error('Uneven number of parameters and values in list.');
        end
        
        tmp = reshape(varargin, [2, nn / 2]);
        for kk = 1 : size(tmp, 2)
            parameters.(tmp{1, kk}) = tmp{2, kk};
        end
    end
end