clc
clear all
%%Nonlinear Test
x = [0:0.4:4,4,4:-0.4:0,0];
y =[zeros([1,11]),.1,.2*ones([1,11]),.1];
ux0 = zeros(size(x));
uy0 = zeros(size(y));
for k = 1:10
    [xo,yo,ux,uy] = sc11(x+ux0,y+uy0);
    ux0 = ux0+ux;
    uy0 = uy+uy0;
end


coe = 10;
plot(x,y,x+coe*ux0,y + coe*uy0);