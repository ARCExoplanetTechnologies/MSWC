function [plane] = create_opticalPlane(N, D, gridCentering)

if strcmp(gridCentering,'vertex')
    centerFactor = 1;
elseif strcmp(gridCentering,'cell')
	centerFactor = 1/2;
end

plane.N = N;
plane.Dx = D; 
plane.Dy = D;
plane.D = D;
plane.dx = plane.Dx/plane.N;
plane.dy = plane.Dy/plane.N;
plane.x = ((1:plane.N) - centerFactor - plane.N/2)/(plane.N)*plane.D;
plane.y = plane.x;
[plane.xx plane.yy] = meshgrid(plane.x, plane.y);
plane.rr = sqrt(plane.xx.^2 + plane.yy.^2);
plane.ttheta = atan2(plane.yy, plane.xx);

end

