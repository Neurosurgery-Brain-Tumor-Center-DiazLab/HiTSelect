function y=obj_fun(A,e,x)
%function y=obj_fun(A,e,x)
%
%IN: A is map from string ids to sparse n-by-n adjacency matrices
%    e is a vector, e(i) is the largest eigenvalue of A(kz{i}) where
%       kz=A.keys
%    x is a vector of length n
%
%OUT: y is the value of an objective function whose extrema defines the
%joint spectrum of the matrices in A

alpha=4;
kz=A.keys;
y=0;%g=0;
for i=1:length(kz)
    v=A(kz{i})*x-e(i)*x;
    y=y+v'*v;
    %g=g+v;
end
y=y+alpha*abs(norm(x)^2-1);
%g=g+alpha*abs(x);