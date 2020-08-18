 function Y=downsample3D(X,a1,a2,a3)
 
X=downsample(X,a1); 
X=permute(X,[2 1 3]);
X=downsample(X,a2);
X=permute(X,[3 1 2]);
X=downsample(X,a3);
Y=permute(X,[3 2 1]);
 end