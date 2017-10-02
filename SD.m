function [ Ds,Ss ] = SD( xi,yi,x1,y1,x2,y2,x3,y3,nu,G )
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
FA = 1-4*nu;
AL = 1-2*nu;


D11 = zeros([1,6]);
D22 = D11;
D12 = D11;
S11 = D11;
S12 = S11;
S22 = S11;


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
    
    %%
    D11(1) = D11(1) + (AL*RDx+2*RDx^3)*gaussw(k)*xja*w1/(DE*RA);
    D11(2) = D11(2) + (2*RDx^2*RDy-AL*RDy)*gaussw(k)*xja*w1/(DE*RA);
    D11(3) = D11(3) + (AL*RDx+2*RDx^3)*gaussw(k)*xja*w2/(DE*RA);
    D11(4) = D11(4) + (2*RDx^2*RDy-AL*RDy)*gaussw(k)*xja*w2/(DE*RA);
    D11(5) = D11(5) + (AL*RDx+2*RDx^3)*gaussw(k)*xja*w3/(DE*RA);
    D11(6) = D11(6) + (2*RDx^2*RDy-AL*RDy)*gaussw(k)*xja*w3/(DE*RA);
    
    D12(1) = D12(1) + (AL*RDy+2*RDx^2*RDy)*w1/(DE*RA)*gaussw(k)*xja;
    D12(2) = D12(2) + (AL*RDx+2*RDy^2*RDx)*w1/(DE*RA)*gaussw(k)*xja;
    D12(3) = D12(3) + (AL*RDy+2*RDx^2*RDy)*w2/(DE*RA)*gaussw(k)*xja;
    D12(4) = D12(4) + (AL*RDx+2*RDy^2*RDx)*w2/(DE*RA)*gaussw(k)*xja;
    D12(5) = D12(5) + (AL*RDy+2*RDx^2*RDy)*w3/(DE*RA)*gaussw(k)*xja;
    D12(6) = D12(6) + (AL*RDx+2*RDy^2*RDx)*w3/(DE*RA)*gaussw(k)*xja;
    
    D22(1) = D22(1) + (2*RDy^2*RDx-AL*RDx)*w1/(DE*RA)*gaussw(k)*xja;
    D22(2) = D22(2) + (AL*RDy + 2*RDy^3)*w1/(DE*RA)*gaussw(k)*xja;
    D22(3) = D22(3) + (2*RDy^2*RDx-AL*RDx)*w2/(DE*RA)*gaussw(k)*xja;
    D22(4) = D22(4) + (AL*RDy + 2*RDy^3)*w2/(DE*RA)*gaussw(k)*xja;
    D22(5) = D22(5) + (2*RDy^2*RDx-AL*RDx)*w3/(DE*RA)*gaussw(k)*xja;
    D22(6) = D22(6) + (AL*RDy + 2*RDy^3)*w3/(DE*RA)*gaussw(k)*xja;

    S11(1) = S11(1) + (2*RDn*(AL*RDx+nu*2*RDx-4*RDx^3)+4*nu*nx*RDx^2+...
        AL*(2*nx*RDx^2+2*nx)-FA*nx)*2*G*w1/(DE*RA^2)*gaussw(k)*xja;
    S11(2) = S11(2) + (2*RDn*(AL*RDy-4*RDx^2*RDy)+4*nu*nx*RDx*RDy+...
        AL*2*ny*RDx^2-FA*ny)*2*G*w1/(DE*RA^2)*gaussw(k)*xja;
    S11(3) = S11(3) + (2*RDn*(AL*RDx+nu*2*RDx-4*RDx^3)+4*nu*nx*RDx^2+...
        AL*(2*nx*RDx^2+2*nx)-FA*nx)*2*G*w2/(DE*RA^2)*gaussw(k)*xja;
    S11(4) = S11(4) + (2*RDn*(AL*RDy-4*RDx^2*RDy)+4*nu*nx*RDx*RDy+...
        AL*2*ny*RDx^2-FA*ny)*2*G*w2/(DE*RA^2)*gaussw(k)*xja;
    S11(5) = S11(5) + (2*RDn*(AL*RDx+nu*2*RDx-4*RDx^3)+4*nu*nx*RDx^2+...
        AL*(2*nx*RDx^2+2*nx)-FA*nx)*2*G*w3/(DE*RA^2)*gaussw(k)*xja;
    S11(6) = S11(6) + (2*RDn*(AL*RDy-4*RDx^2*RDy)+4*nu*nx*RDx*RDy+...
        AL*2*ny*RDx^2-FA*ny)*2*G*w3/(DE*RA^2)*gaussw(k)*xja;
    
    S12(1) = S12(1) + (2*RDn*(nu*RDy-4*RDx^2*RDy)+2*nu*(nx*RDy*RDx+ny*RDx^2)...
        +AL*(2*nx*RDx*RDy+ny))*2*G*w1/(DE*RA^2)*gaussw(k)*xja;
    S12(2) = S12(2) + (2*RDn*(nu*RDx-4*RDy^2*RDx)+2*nu*(nx*RDy^2+ny*RDx*RDy)...
        +AL*(2*ny*RDx*RDy+nx))*2*G*w1/(DE*RA^2)*gaussw(k)*xja;
    S12(3) = S12(3) + (2*RDn*(nu*RDy-4*RDx^2*RDy)+2*nu*(nx*RDy*RDx+ny*RDx^2)...
        +AL*(2*nx*RDx*RDy+ny))*2*G*w2/(DE*RA^2)*gaussw(k)*xja;
    S12(4) = S12(4) + (2*RDn*(nu*RDx-4*RDy^2*RDx)+2*nu*(nx*RDy^2+ny*RDx*RDy)...
        +AL*(2*ny*RDx*RDy+nx))*2*G*w2/(DE*RA^2)*gaussw(k)*xja;
    S12(5) = S12(5) + (2*RDn*(nu*RDy-4*RDx^2*RDy)+2*nu*(nx*RDy*RDx+ny*RDx^2)...
        +AL*(2*nx*RDx*RDy+ny))*2*G*w3/(DE*RA^2)*gaussw(k)*xja;
    S12(6) = S12(6) + (2*RDn*(nu*RDx-4*RDy^2*RDx)+2*nu*(nx*RDy^2+ny*RDx*RDy)...
        +AL*(2*ny*RDx*RDy+nx))*2*G*w3/(DE*RA^2)*gaussw(k)*xja;
    
    S22(1) = S22(1) + (2*RDn*(AL*RDx-4*RDx*RDy^2)+4*nu*ny*RDx*RDy+...
        AL*2*nx*RDy^2-FA*nx)*2*G*w1/(DE*RA^2)*gaussw(k)*xja;
    S22(2) = S22(2) + (2*RDn*(AL*RDy+2*nu*RDy-4*RDy^3)+4*nu*ny*RDy^2+...
        AL*(2*ny*RDy^2+2*ny)-FA*ny)*2*G*w1/(DE*RA^2)*gaussw(k)*xja;
    S22(3) = S22(3) + (2*RDn*(AL*RDx-4*RDx*RDy^2)+4*nu*ny*RDx*RDy+...
        AL*2*nx*RDy^2-FA*nx)*2*G*w2/(DE*RA^2)*gaussw(k)*xja;
    S22(4) = S22(4) + (2*RDn*(AL*RDy+2*nu*RDy-4*RDy^3)+4*nu*ny*RDy^2+...
        AL*(2*ny*RDy^2+2*ny)-FA*ny)*2*G*w2/(DE*RA^2)*gaussw(k)*xja;
    S22(5) = S22(5) + (2*RDn*(AL*RDx-4*RDx*RDy^2)+4*nu*ny*RDx*RDy+...
        AL*2*nx*RDy^2-FA*nx)*2*G*w3/(DE*RA^2)*gaussw(k)*xja;
    S22(6) = S22(6) + (2*RDn*(AL*RDy+2*nu*RDy-4*RDy^3)+4*nu*ny*RDy^2+...
        AL*(2*ny*RDy^2+2*ny)-FA*ny)*2*G*w3/(DE*RA^2)*gaussw(k)*xja;
end
Ds = [D11;D12;D22];
Ss = [S11;S12;S22];

end

