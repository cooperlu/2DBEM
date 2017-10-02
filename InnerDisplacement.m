function [ uix,uiy ] = InnerDisplacement( x,y,xi,yi,ux,uy,n,nu,mu,f,u )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
ne = n/2;
ll=1;
[in,on] = inpolygon(xi,yi,x,y);
if on == 0    
    Hi = zeros([2,4*ne]);
    Gi = zeros([2,6*ne]);
    
    for i = 1:2:n-1
        [GWs,HWs] = HGndm(xi,yi,x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
        
        Gi(2-1:2,3*i-2:3*i+3) = Gi(2-1:2,3*i-2:3*i+3) + GWs;
        if i == n-1
            Hi(2*ll-1:2*ll,2*i-1:2*i+2) = Hi(2*ll-1:2*ll,2*i-1:2*i+2) + HWs(1:2,1:4);
            Hi(2*ll-1:2*ll,1:2) = Hi(2*ll-1:2*ll,1:2) + HWs(1:2,5:6);
        else
            Hi(2*ll-1:2*ll,2*i-1:2*i+4) = Hi(2*ll-1:2*ll,2*i-1:2*i+4) + HWs;
        end
    end
    
    ui = (Gi*f - Hi * u)';
    uix = ui(1);
    uiy = ui(2);
    return
else
    xo = x(1:end-1);
    yo = y(1:end-1);
    if sum((xi == xo).*(yi ==yo)) > 0
        uix = ((xi == xo).*(yi ==yo)) * ux';
        uiy = ((xi == xo).*(yi ==yo)) * uy';
        return
    else
        uxm = ux;
        uxm(end+1) = ux(1);
        uym = uy;
        uym(end+1) = uy(1);
        for i = 1:ne
            x0 = x(i*2-1);
            y0 = y(i*2-1);
            x1 = x(i*2);
            y1 = y(i*2);
            x2 = x(i*2+1);
            y2 = y(i*2+1);
            xp0 = xi-x0;
            yp0 = yi-y0;
            xp1 = xi-x1;
            yp1 = yi-y1;
            xp2 = xi-x2;
            yp2 = yi-y2;
            if (xp0*yp2 == xp2*yp0)&&((xp0*xp2<0)||(yp0*yp2<0))
                if(x2-x1) ==0
                    z = yp1/(y2-y1);
                else
                    z = xp1/(x2-x1);
                end
                w1 = z*(z-1)*0.5;
                w2 = 1-z^2;
                w3 = z*(z+1)*0.5;
                            ux0 = uxm(i*2-1);
            uy0 = uym(i*2-1);
            ux1 = uxm(i*2);
            uy1 = uym(i*2);
            ux2 = uxm(i*2+1);
            uy2 = uym(i*2+1);
            uix = w1*ux0+w2*ux1+w3*ux2;
            uiy = w1*uy0+w2*uy1+w3*uy2;
            return
            end
        end
    end
end
end

