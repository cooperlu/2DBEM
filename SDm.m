function [ Ds,Ss ] = SD( xi,yi,x1,y1,x2,y2,x3,y3,nu,G )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Gauss Quadrature Value
Gauss =[1	0.0486909570091397	-0.0243502926634244;
2	0.0486909570091397	0.0243502926634244;
3	0.0485754674415034	-0.0729931217877990;
4	0.0485754674415034	0.0729931217877990;
5	0.0483447622348030	-0.1214628192961206;
6	0.0483447622348030	0.1214628192961206;
7	0.0479993885964583	-0.1696444204239928;
8	0.0479993885964583	0.1696444204239928;
9	0.0475401657148303	-0.2174236437400071
10	0.0475401657148303	0.2174236437400071;
11	0.0469681828162100	-0.2646871622087674;
12	0.0469681828162100	0.2646871622087674;
13	0.0462847965813144	-0.3113228719902110;
14	0.0462847965813144	0.3113228719902110;
15	0.0454916279274181	-0.3572201583376681;
16	0.0454916279274181	0.3572201583376681;
17	0.0445905581637566	-0.4022701579639916;
18	0.0445905581637566	0.4022701579639916;
19	0.0435837245293235	-0.4463660172534641;
20	0.0435837245293235	0.4463660172534641;
21	0.0424735151236536	-0.4894031457070530;
22	0.0424735151236536	0.4894031457070530;
23	0.0412625632426235	-0.5312794640198946;
24	0.0412625632426235	0.5312794640198946;
25	0.0399537411327203	-0.5718956462026340;
26	0.0399537411327203	0.5718956462026340;
27	0.0385501531786156	-0.6111553551723933;
28	0.0385501531786156	0.6111553551723933;
29	0.0370551285402400	-0.6489654712546573;
30	0.0370551285402400	0.6489654712546573;
31	0.0354722132568824	-0.6852363130542333;
32	0.0354722132568824	0.6852363130542333;
33	0.0338051618371416	-0.7198818501716109;
34	0.0338051618371416	0.7198818501716109;
35	0.0320579283548516	-0.7528199072605319;
36	0.0320579283548516	0.7528199072605319;
37	0.0302346570724025	-0.7839723589433414;
38	0.0302346570724025	0.7839723589433414;
39	0.0283396726142595	-0.8132653151227975;
40	0.0283396726142595	0.8132653151227975;
41	0.0263774697150547	-0.8406292962525803;
42	0.0263774697150547	0.8406292962525803;
43	0.0243527025687109	-0.8659993981540928;
44	0.0243527025687109	0.8659993981540928;
45	0.0222701738083833	-0.8893154459951141;
46	0.0222701738083833	0.8893154459951141;
47	0.0201348231535302	-0.9105221370785028;
48	0.0201348231535302	0.9105221370785028;
49	0.0179517157756973	-0.9295691721319396;
50	0.0179517157756973	0.9295691721319396;
51	0.0157260304760247	-0.9464113748584028;
52	0.0157260304760247	0.9464113748584028;
53	0.0134630478967186	-0.9610087996520538;
54	0.0134630478967186	0.9610087996520538;
55	0.0111681394601311	-0.9733268277899110;
56	0.0111681394601311	0.9733268277899110;
57	0.0088467598263639	-0.9833362538846260;
58	0.0088467598263639	0.9833362538846260;
59	0.0065044579689784	-0.9910133714767443;
60	0.0065044579689784	0.9910133714767443;
61	0.0041470332605625	-0.9963401167719553;
62	0.0041470332605625	0.9963401167719553;
63	0.0017832807216964	-0.9993050417357722;
64	0.0017832807216964	0.9993050417357722];
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


for k = 1:64
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

