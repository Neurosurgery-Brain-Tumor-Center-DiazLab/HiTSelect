function z=grad_fun(A,e,x)
%function z=grad_fun(A,e,x)
%
%IN: A is map from string ids to sparse n-by-n adjacency matrices
%    e is a vector, e(i) is the largest eigenvalue of A(kz{i}) where
%       kz=A.keys
%    x is a vector of length n
%
%OUT: y is the value of the gradient of obj_fun.m, the objective function whose extrema defines the
%joint spectrum of the matrices in A

kz=A.keys;
y=0;
for i=1:length(kz)
    v=A(kz{i})*x-e(i)*x;
    y=y+v;
end
y=2*(y+abs(x));