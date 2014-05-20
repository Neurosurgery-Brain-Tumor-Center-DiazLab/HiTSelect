function [c,ceq] = mycon(x)
%function [c,ceq] = mycon(x)
%
%IN: x is a vector
%
%OUT: c==[] and ceq==norm(x)-1

c=[];ceq=abs(norm(x)-1);