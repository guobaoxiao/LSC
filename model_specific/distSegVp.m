function [ d ] = distSegVp( x, V )
%DISTSEGVP compute the distance between a point and a vanishing point
centroid = [x(1,:)+x(2,:); x(3,:)+x(4,:)]./2;
l=cross([centroid;1], V);
d = distPointLine([x(1,:),x(3,:)]', l);

end

