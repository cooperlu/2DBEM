function [ Gs,Hs ] = HGnd( xi,yi,x1,y1,x2,y2,x3,y3,nu,G )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Gauss Quadrature Value
Gauss =[1	0.2955242247147529	-0.1488743389816312;
2	0.2955242247147529	0.1488743389816312;
3	0.2692667193099963	-0.4333953941292472;
4	0.2692667193099963	0.4333953941292472;
5	0.2190863625159820	-0.6794095682990244;
6	0.2190863625159820	0.6794095682990244;
7	0.1494513491505806	-0.8650633666889845;
8	0.1494513491505806	0.8650633666889845;
9	0.0666713443086881	-0.9739065285171717;
10	0.0666713443086881	0.9739065285171717;];
gaussw = Gauss(:,2);
gaussp = Gauss(:,3);

%% Geometric value

A = x3-2*x2+x1;
B = (x3-x1)/2;
C = y3-2*y2+y1;
D = (y3-y1)/2;
DE = 4*pi*(1-nu);


% normal = [direction(2),-direction(1)];
% lengthj = sqrt(direction(1)^2 + direction(2)^2);
% nx = normal(1)/lengthj;
% ny = normal(2)/lengthj;

Gs = zeros([2,6]);
Hs = zeros([2,6]);

for k = 1:10
    w1 = gaussp(k)*(gaussp(k)-1)*0.5;
    w2 = 1-gaussp(k)^2;
    w3 = gaussp(k)*(gaussp(k)+1)*0.5;
    W = [w1,0,w2,0,w3,0;
        0,w1,0,w2,0,w3];
    XCO = x1*w1 + x2*w2 + x3*w3;
    YCO = y1*w1 + y2*w2 + y3*w3;
    xja = sqrt((gaussp(k)*A + B)^2 + (gaussp(k)*C + D)^2);
    nx = (gaussp(k)*C + D)/xja;
    ny = -(gaussp(k)*A + B)/xja;
    
    RA = sqrt((xi-XCO)^2 + (yi-YCO)^2);
    RDx = (XCO - xi)/RA;
    RDy = (YCO - yi)/RA;
    RDn = RDx*nx + RDy*ny;
    
    
    G11 = ((3-4*nu)*log(1/RA) + RDx^2)*gaussw(k)*xja ...
        /(2*DE*G);
    G12 = RDx*RDy*gaussw(k)*xja/(2*DE*G);
    G22 = ((3-4*nu)*log(1/RA) + RDy^2)*gaussw(k)*xja ...
        /(2*DE*G);
    H11 = - RDn*((1-2*nu) + 2*RDx^2)/(RA*DE)*gaussw(k)...
        *xja;
    H12 = - (RDn*2*RDx*RDy + (1-2*nu)*(nx*RDy - ny*RDx))...
        *gaussw(k)*xja/(RA*DE);
    H21 = - (RDn*2*RDx*RDy - (1-2*nu)*(nx*RDy - ny*RDx))...
        *gaussw(k)*xja/(RA*DE);
    H22 = - RDn*((1-2*nu) + 2*RDy^2)/(RA*DE)*gaussw(k)...
        *xja;
    Gs = Gs + [G11,G12;
               G12,G22]*W;
    Hs = Hs + [H11,H12;
               H21,H22]*W;
    
end


end

