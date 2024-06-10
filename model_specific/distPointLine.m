function [ d ] = distPointLine( P,l )
%DISTPOINTLINE compute the euclidean distance between a point P and a line L in R^2.
%   Detailed explanation goes here

d=abs(l(1).*P(1,:)+l(2).*P(2,:)+l(3))/sqrt(l(1)^2+l(2)^2);

end

