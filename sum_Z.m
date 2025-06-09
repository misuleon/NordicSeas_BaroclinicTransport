function ans = sum_Z(z,prop)
% routine sums property in z
% z must be column vector
[m,n] = size(prop);
dz = diff(z);
DZ = repmat(dz,1,n);
prop_2 = (prop(2:m,:) + prop(1:m-1,:))/2;
ans = sum(DZ.*prop_2,1);
return
