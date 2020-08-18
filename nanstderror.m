function a0=nanstderror(a,n)
if nargin == 1
a0=nanstd(a)/sqrt(sum(~isnan(a))-1);

else
a0=nanstd(a,n)./sqrt(sum(~isnan(a),n)-1);
end