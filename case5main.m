% function [xo,yo,ux,uy] = sc11(x,y)
clc
clear all

%% Parameters
mu = 0.8e5;
nu = 0.25;
n = 30;
nn = 2*n; %number of unknows
ne = n/2; %number of elements



%% Geometry
x = [-5:0.25:-1,-sqrt(2)/2,zeros([1,4]),0:-1.25:-5,-5*ones([1,3])];
y = [zeros([1,17]),sqrt(2)/2,1:4:17,17*ones([1,3]),17:-4.25:4.25];
xo = x;
yo = y;
x(1,end+1) = x(1);
y(1,end+1) = y(1);


%% BC

f = zeros([3*n,1]);
f(68:2:78) = 100;

bct = ones([3*n,1]);
dc = [2:2:48,55:2:65];
bct(dc) = 0;


%% BUILD H & G

[ Hma,Gma ] = HG( x,y,nu,mu,n,ne );


%% Apply BC
H = Hma;
G = Gma;

[ H,G ] = applybc( mu,ne,bct,H,G );
                  
%% SOLVE

u = H\(G*f);

%% Organize


u = u;
f = f;

for i = 1:ne
    for j = 1:6 % 6 unknowns each element
        if bct(6*i-6+j) == 0 % 0 is displacement, 1 is traction
            if (i == ne) && j > 4 % 1st code in the last element
                if bct(j-4) > 0
                    ch = u(j-4);
                    u(j-4) =  f(6*i-6+j);
                    f(6*i-6+j) = ch * mu;
                else
                    f(6*i-6+j) = f(j-4);
                end
            else
                if (i == 1) || (j>2) || (bct(6*i-8+j) ==1)
                    ch = u(4*i-4+j);
                    u(4*i-4+j) =  f(6*i-6+j);
                    f(6*i-6+j) = ch * mu;
                else
                    f(6*i-6+j) =f(6*i-8+j);
                end
            end
        end
    end
end
%% Calculate internal points

% 
% for i = 1:19
%     for j = 1:79
%         xi((i-1)*79+j) = j*0.05;
%         yi((i-1)*79+j) = i*0.01;
%     end
% end
xi = -1;
yi = 0.1;

[~,xs] = size(xi);

Hi = zeros([2*xs,4*ne]);
Gi = zeros([2*xs,6*ne]);
Di = zeros([3*xs,6*ne]);
Si = zeros([3*xs,4*ne]);
for ll = 1:xs
for i = 1:2:n-1
[GWs,HWs] = HGnd(xi(ll),yi(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
        [ Ds,Ss ] = SDm( xi(ll),yi(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu );
        Gi(2*ll-1:2*ll,3*i-2:3*i+3) = Gi(2*ll-1:2*ll,3*i-2:3*i+3) + GWs;
Di(3*ll-2:3*ll,3*i-2:3*i+3) = Di(3*ll-2:3*ll,3*i-2:3*i+3) + Ds;
        if i == n-1
            Hi(2*ll-1:2*ll,2*i-1:2*i+2) = Hi(2*ll-1:2*ll,2*i-1:2*i+2) + HWs(1:2,1:4);
            Hi(2*ll-1:2*ll,1:2) = Hi(2*ll-1:2*ll,1:2) + HWs(1:2,5:6);
            Si(3*ll-2:3*ll,2*i-1:2*i+2) = Si(3*ll-2:3*ll,2*i-1:2*i+2) + Ss(1:3,1:4);
            Si(3*ll-2:3*ll,1:2) = Si(3*ll-2:3*ll,1:2) + Ss(1:3,5:6);
        else
            Hi(2*ll-1:2*ll,2*i-1:2*i+4) =Hi(2*ll-1:2*ll,2*i-1:2*i+4) + HWs;
            Si(3*ll-2:3*ll,2*i-1:2*i+4) =Si(3*ll-2:3*ll,2*i-1:2*i+4) + Ss;
        end
end
end
ui = Gi*f - Hi * u;
sigmai = Di*f-Si*u;
for i = 1:xs
    uix = ui((i-1)*2+1);
    uiy = ui(i*2);
    uia(i) = sqrt(uix^2 + uiy^2);
end

% for i = 1:19
%     for j = 1:79
%         xim(i,j) = xi((i-1)*79+j);
%         yim(i,j) = yi((i-1)*79+j);
%         uim(i,j) = uia((i-1)*79+j);
%     end
% end

% surf(xim,yim,uim)
ux = 0;
uy = 0;
for i = 1:n
    ux(i) = u(2*i-1);
    uy(i) = u(2*i);
end
for i = 1:n
    fx(i) = f(2*i-1);
    fy(i) = f(2*i);
end
ux(end+1) = ux(1);
uy(end+1) = uy(1);
coe = 10;
plot(x,y,x+ux*coe,y+uy*coe);