function [] = display_correspondence(img1,img2, X,G )
%DISPLAY_CORRESPONDENCE Summary of this function goes here
%   Detailed explanation goes here
%%
[h]=size(img1,1);
cmap=colormap('Lines');
im=[img1;img2];
N=size(X,2);
Y=X;
figure
imshow(im); hold all;
Y(5,:)=Y(5,:)+h;
gscatter(Y(1,:),Y(2,:),G);
gscatter(Y(4,:),Y(5,:),G); legend off;
line([0,0], [h,h])
G=G+1;
for i=1:N
    line(Y([1,4],i),Y([2,5],i),'Color',cmap(G(i),:))
end
end

