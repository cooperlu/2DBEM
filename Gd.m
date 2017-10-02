function [GWs] = Gd(xg1,yg1,xg2,yg2,xg3,yg3,nu,G,P )
% Calculate G for local points
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
gaussw = (Gauss(:,2));
gaussp = (Gauss(:,3));

SQp = [0.0090426309,0.0539712662,0.1353118246,0.2470524162,0.3802125396,0.5237923179,...
    0.6657752055,0.7941904160,0.8981610912,0.9688479887];
SQw = [0.1209551319,0.1863635425,0.1956608732,0.1735771421,0.1356956729,0.0936467585,...
    0.0557877273,0.0271598109,0.0095151826,0.0016381576];

%% Geometric value
switch P
    case 1
        x3 = xg3 - xg1;
        y3 = yg3 - yg1;
        x2 = xg2 - xg1;
        y2 = yg2 - yg1;
        a1 = (x3-2*x2)*0.5;
        b1 = x2;
        a2 = (y3-2*y2)*0.5;
        b2 = y2;
    case 2
        x3 = xg3 - xg2;
        y3 = yg3 - yg2;
        x1 = xg1 - xg2;
        y1 = yg1 - yg2;
        a1 = x1 + x3;
        b1 = x3 - x1;
        a2 = y1 + y3;
        b2 = y3 - y1;
    case 3
        x2 = xg2 - xg3;
        y2 = yg2 - yg3;
        x1 = xg1 - xg3;
        y1 = yg1 - yg3;
        a1 = (x1-2*x2)*0.5;
        b1 = -x2;
        a2 = (y1-2*y2)*0.5;
        b2 = -y2;
end
GWs = zeros([2,6]);


A = a1^2 + a2^2;
B = 2*(a1*b1 +a2*b2);
C = b1^2 + b2^2;
cont1 = (3-4*nu)/(8*pi*G*(1-nu));
cont2 = 1/(8*pi*G*(1-nu));

%% Loop
for k = 1:10
    T1 = ((a1*gaussp(k) + b1)^2)/((a2*gaussp(k) + b2)^2 +(a1*gaussp(k) + b1)^2);
    T2 = ((a2*gaussp(k) + b2)^2)/((a2*gaussp(k) + b2)^2 +(a1*gaussp(k) + b1)^2);
    T3 = ((a2*gaussp(k) + b2)*(a1*gaussp(k) + b1))/((a2*gaussp(k) + b2)^2 +(a1*gaussp(k) + b1)^2);

    %Shape Function for Integration

    f1 = gaussp(k)*(gaussp(k)-1)*0.5;
    f2 = 1-gaussp(k)^2;
    f3 = gaussp(k)*(gaussp(k)+1)*0.5;
    fl1 = (SQp(k)-1)*(2*SQp(k)-1);
    fl2 = 4* SQp(k)*(1-SQp(k));
    fl3 = (SQp(k))*(2*SQp(k)-1);
    fln1 = SQp(k)*(SQp(k)-1)*0.5;
    fln2 = 1 - SQp(k)^2;
    fln3 = SQp(k)*(SQp(k) +1)*0.5;
    
    % Compute GW
    switch P
        case 1
            xja1 = sqrt((4*a1*SQp(k)-2*a1+0.5*x3)^2+(4*a2*SQp(k)-2*a2 +...
                0.5*y3)^2)*2;
            xja2 = sqrt((a1*gaussp(k)*2+0.5*x3)^2+(a2*gaussp(k)*2+0.5*y3)^2);
            xlo = -log(2*sqrt((gaussp(k)*a1+b1)^2 + (gaussp(k)*a2 + b2)^2));
            s1 = cont1*(fl1*xja1*SQw(k) + f1 * xja2 * xlo * gaussw(k));
            s2 = cont1*(fl2*xja1*SQw(k) + f2 * xja2 * xlo * gaussw(k));
            s3 = cont1*(fl3*xja1*SQw(k) + f3 * xja2 * xlo * gaussw(k));
        case 2
            xja1 = sqrt((0.5*b1 - a1*SQp(k))^2 + (0.5*b2 - a2*SQp(k))^2);
            xja11= sqrt((0.5*b1 + a1*SQp(k))^2 + (0.5*b2 + a2*SQp(k))^2);
            xja2 = sqrt((0.5*b1 + a1*gaussp(k))^2 + (0.5*b2 + a2*gaussp(k))^2);
            xlo = -0.5*log((gaussp(k)*a1*0.5 + b1*0.5)^2 +(gaussp(k)*a2*0.5 + b2*0.5)^2);
            s1 = cont1*((fln3*xja1 + fln1*xja11)*SQw(k) + f1*xja2*xlo*gaussw(k));
            s2 = cont1*(fln2*(xja1+xja11)*SQw(k) + f2*xja2*xlo*gaussw(k));
            s3 = cont1*((fln1*xja1 + fln3*xja11)*SQw(k) + f3*xja2*xlo*gaussw(k));
        case 3
            xja1 = sqrt((2*a1 -4*a1*SQp(k)-0.5*x1)^2+(-4*a2*SQp(k)+2*a2 -...
                0.5*y1)^2)*2;
            xja2 = sqrt((a1*gaussp(k)*2-0.5*x1)^2+(a2*gaussp(k)*2-0.5*y1)^2);
            xlo = -log(2*sqrt((gaussp(k)*a1+b1)^2 + (gaussp(k)*a2 + b2)^2));
            s1 = cont1*(fl3*xja1*SQw(k) + f1 * xja2 * xlo * gaussw(k));
            s2 = cont1*(fl2*xja1*SQw(k) + f2 * xja2 * xlo * gaussw(k));
            s3 = cont1*(fl1*xja1*SQw(k) + f3 * xja2 * xlo * gaussw(k));
    end

    GWs(1,1) = GWs(1,1) + s1 + cont2*f1*T1*xja2*gaussw(k);
    GWs(1,3) = GWs(1,3) + s2 + cont2*f2*T1*xja2*gaussw(k);
    GWs(1,5) = GWs(1,5) + s3 + cont2*f3*T1*xja2*gaussw(k);
    GWs(1,2) = GWs(1,2) + cont2*f1*T3*xja2*gaussw(k);
    GWs(1,4) = GWs(1,4) + cont2*f2*T3*xja2*gaussw(k);
    GWs(1,6) = GWs(1,6) + cont2*f3*T3*xja2*gaussw(k);
    GWs(2,1) = GWs(1,2);
    GWs(2,3) = GWs(1,4);
    GWs(2,5) = GWs(1,6);
    GWs(2,2) = GWs(2,2) + s1 + cont2*f1*T2*xja2*gaussw(k);
    GWs(2,4) = GWs(2,4) + s2 + cont2*f2*T2*xja2*gaussw(k);
    GWs(2,6) = GWs(2,6) + s3 + cont2*f3*T2*xja2*gaussw(k);
    
end


end

