# honeycomb
> Tools for hexagonal binning (honeycomb plot) and visualisation

## Installation
Put all files or at least `honeycomb.m` somewhere on your MATLAB path.

## Documentation
Run `help honeycomb` or `doc honeycomb` for instructions.

## Example
```matlab
x = randn(100);
y = rand(size(x));
figure
honeycomb(x, y)
colorbar
title 'Honeycomb plot of uniform vs. normal random data'
```

Result:

![honeycomb](https://cloud.githubusercontent.com/assets/19374736/24495979/3f70f6f6-1537-11e7-8203-523ae248bffa.png)
