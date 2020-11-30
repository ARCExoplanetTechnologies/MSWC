function [outx,outy] = rotcoords2d(theta,inx,iny)
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    
    coordsin = [inx; iny]; 
    coordsout = R*coordsin;
    
    outx = coordsout(1);
    outy = coordsout(2);
end

