% Unit tests for hexscatter-bean

clear
close all

%% hexscatter(x, y)
% Test happy path
x = randn(100, 100);
y = rand(size(x));
figure
hexscatter(x, y)

%TODO
% Assert errors upon any invalid x and or y
% Assert empty axes upon x = [] and y = [], compare: scatter([], [])

%% hexscatter(x, y, 'XLim', XLim)
XLim = [-10, 10];
figure
hexscatter(x, y, 'XLim', XLim)

%% hexscatter(x, y, 'YLim', YLim)
YLim = [-0.1, 1.1];
figure
hexscatter(x, y, 'YLim', YLim)

%% hexscatter(x, y, 'Resolution', Resolution)
Resolution = 2;
figure
try
    hexscatter(x, y, 'Resolution', Resolution)
catch
    disp FIXME! %FIXME
end

%% hexscatter(x, y, 'DrawEdges', DrawEdges)
DrawEdges = true;
figure
hexscatter(x, y, 'DrawEdges', DrawEdges)

%% hexscatter(x, y, 'ShowZeros', ShowZeros)
ShowZeros = true;
figure
hexscatter(x, y, 'ShowZeros', ShowZeros)

%% h = hexscatter(x, y)
figure
h = hexscatter(x, y);
assert(isa(h, 'matlab.graphics.primitive.Patch'))