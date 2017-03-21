function h = hexscatter(xdata, ydata, varargin)
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
%   Example: 'xlim', [0, 1], 'ylim', [1, 10]
%
%   'xlim': x-axis limits
%   [min(x(:)), max(x(:))] (default) | two-element numeric vector
%   The x-axis limits specified as a two-element numeric vector. The first
%   and second element correspond to the left and right x-axis limit
%   respectively. The second element must be greater than the first.
%
%   'ylim': y-axis limits
%   [min(y(:)), max(y(:))] (default) | two-element numeric vector
%   The y-axis limits specified as a two-element numeric vector. The first
%   and second element correspond to the lower and upper y-axis limit
%   respectively. The second element must be greater than the first.
%
%   'res': resolution
%   [100, 100] (default) | numeric scalar | two-element numeric vector
%   The hexagonal grid resolution specified either as a scalar, or a
%   two-element vector. The resolution is the number of hexagons in
%   horizontal and vertical direction. If the specified value is a
%   two-element vector, the first element corresponds to the vertical
%   resolution and the second to the horizontal resolution. If the
%   specified value is scalar, then its value is used for both directions.
%
%   'drawEdges': draw edges around hexagons
%   false (default) | true
%   Draw edges around the hexagons if the specified value is true.
%
%   'showZeros': show zero count hexagons
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


%% Convert to vectors
xdata = xdata(:);
ydata = ydata(:);

nix = isnan(xdata) | isnan(ydata) | isinf(xdata) | isinf(ydata);
xdata = xdata(~nix);
ydata = ydata(~nix);

params = default_param( varargin, ...
    'xlim', [min(xdata(:)) max(xdata(:))], ...
    'ylim', [min(ydata(:)) max(ydata(:))], ...
    'res', [100 100], ...
    'drawEdges', false, ...
    'showZeros', false);

if params.drawedges
    ec = 'flat';
else
    ec = 'none';
end

if length(params.res) == 1
    params.res = params.res * [1, 1];
end

%% Determine grid
xl = params.xlim;
yl = params.ylim;

xbins = linspace(xl(1), xl(2), params.res(2));
ybins = linspace(yl(1), yl(2), params.res(1));
dy = diff(ybins([1, 2])) * 0.5;

[X, Y] = meshgrid(xbins, ybins);
[n, m] = size(X);
Y(:, 1 : fix(end / 2) * 2) = ...
    Y(:, 1 : fix(end / 2) * 2) + repmat([0, dy], [n, fix(m / 2)]);

%% Map points to boxes
% Which pair of columns?
dx = diff(xbins([1, 2]));
foox = floor((xdata - xbins(1)) ./ dx) + 1;
foox(foox > length(xbins)) = length(xbins);
foox(foox < 1) = 1;

% Which pair of rows?
% Use the first row, which starts without an offset, as the standard
fooy = floor((ydata - ybins(1)) ./ diff(ybins([1, 2]))) + 1;
fooy(fooy > length(ybins)) = length(ybins);
fooy(fooy < 1) = 1;

% Which orientation
orientation = mod(foox, 2) == 1;

% Map points to boxes
xdata = xdata - xbins(foox).';
ydata = ydata - ybins(fooy).';

% Which layer
layer = ydata > dy;

% Convert to block B format
toflip = layer == orientation;
xdata(toflip) = dx - xdata(toflip);

ydata(layer == 1) = ydata(layer == 1) - dy;

% Find closest corner
topright = xdata.^2 + ydata.^2 > (xdata - dx).^2 + (ydata - dy).^2;
clear xdata ydata

%% Map corners back to bins
% Which x bin?
foox = foox + ~(orientation == (layer == topright));
foox(foox > length(xbins)) = length(xbins);

% Which y bin?
fooy = fooy + (layer & topright);
fooy(fooy > length(ybins)) = length(ybins);

foox = sub2ind(size(X), fooy, foox);
% Replaced foox to conserve memory

%% Determine counts
counts = histc(foox, 1 : numel(X)).';

newplot;
xscale = diff(xbins([1, 2])) * 2 / 3;
yscale = diff(ybins([1, 2])) * 2 / 3;
theta = 0 : 60 : 360;
x = bsxfun(@plus, X(:), cosd(theta) * xscale).';
y = bsxfun(@plus, Y(:), sind(theta) * yscale).';

if params.showzeros
    h = patch(x, y, counts, 'edgeColor', ec);
else
    jj = counts > 0;
    h = patch(x(:, jj), y(:, jj), counts(jj), 'edgeColor', ec);
end

xlim(params.xlim);
ylim(params.ylim);

if nargout == 0
    clear h;
end

%% Function: default_param
% Gordon Bean, March 2012
% Copied from https://github.com/brazilbean/bean-matlab-toolkit
    function params = default_param(params, varargin)
        if (iscell(params))
            params = get_params(params{:});
        end
        defaults = get_params(varargin{:});
        
        for f = fieldnames(defaults).'
            field = f{:};
            if (~isfield(params, lower(field)))
                params.(lower(field)) = defaults.(field);
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
    function params = get_params(varargin)
        params = struct;
        
        nn = length(varargin);
        if (mod(nn, 2) ~= 0)
            error('Uneven number of parameters and values in list.');
        end
        
        tmp = reshape(varargin, [2, nn / 2]);
        for kk = 1 : size(tmp, 2)
            params.(lower(tmp{1, kk})) = tmp{2, kk};
        end
    end
end