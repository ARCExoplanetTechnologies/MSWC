function [Out,vect]=AZAV(I,ratio,OWA)

L=length(I);
u=(-L/2+0.5:L/2-0.5)/L*2*OWA;
[U,V]=meshgrid(u);
[T,R]=cart2pol(U,V);

%AZAV returns Out, the azymutal average of the 2D array I, and vect, which
% is the new radial vector.
% - ratio is the ratio of the number of points along a radius in the original 
% domain (half of the length of I) and the number of points in vect. This 
% ratio should be small enough to get enough points along a radius.
% - OWA is the physical size of the radius (typically the outer working angle)

M=floor(length(I)/2/ratio);

Out=zeros(M,1);

vect=(0.5:M-0.5)/M*OWA;

for k=1:M
    
    REG=(R<=(vect(k)+ OWA/(2*M))).*(R>(vect(k)- OWA/(2*M))); 
    Out(k)=mean(I(REG==1));
    
end