function [ Y ] = inversegampdf_tune( X,A,B )
%inversegampdf Inverse gamma probability density function.
%   Y = inversegampdf(X,A,B) returns the inverse gamma probability density
%   function with shape and scale parameters A and B, respectively, at the
%   values in X. The size of Y is the common size of the input arguments. A
%   scalar input functions is a constant matrix of the same size as the
%   other inputs.It is modified to have vector as beta. 

if size(X,1) ~= size(B,1)
    X = X';
end

Y = B.^A./gamma(A).*X.^(-A-1).*exp(-B./X);