function [refPosw,wdA] = allrefpoints(w,l,h,r)
%[X,Y,Z] = ndgrid(0.1:0.1:4.9);
%plot3( X(:), Y(:), Z(:))
wdA = r^2;
numpoints = 10000;
 x = rand(numpoints,1) * w;
 y = rand(numpoints,1) * l;
 z = rand(numpoints,1) * h;
refPosw = [x, y, z];
scatter3(x, y, z)
axis([0 w 0 l 0 h])
            axis square
            grid on
            xlabel('x')
            ylabel('y')
            zlabel('z')
            title('position of the reflection points')
            drawnow