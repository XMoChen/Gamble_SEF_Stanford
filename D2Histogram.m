function H=D2Histogram(X0,xbins,ybins)
X = X0(:,1);%randn(2500,1);
Y = X0(:,2);%randn(2500,1)*2;

%# bin centers (integers)
% xbins = floor(min(X)):1:ceil(max(X));
% ybins = floor(min(Y)):1:ceil(max(Y));
xNumBins = numel(xbins); yNumBins = numel(ybins);

%# map X/Y values to bin indices
Xi = round( interp1(xbins, 1:xNumBins, X, 'linear', 'extrap') );
Yi = round( interp1(ybins, 1:yNumBins, Y, 'linear', 'extrap') );

%# limit indices to the range [1,numBins]
Xi = max( min(Xi,xNumBins), 1);
Yi = max( min(Yi,yNumBins), 1);

%# count number of elements in each bin
H = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);