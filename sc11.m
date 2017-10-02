function [xo,yo,ux,uy] = sc11(x,y)


%% Parameters
mu = 0.945e5;
nu = 0.2;
n = 24;
nn = 2*n; %number of unknows
ne = n/2; %number of elements



%% Geometry
xo = x;
yo = y;
x(1,end+1) = x(1);
y(1,end+1) = y(1);


%% BC

f = zeros([3*n,1]);
f(31:2:35) = 100;

bct = ones([3*n,1]);
dc = [67:72];
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
% 
% 
% [~,xs] = size(xi);
% 
% Hi = zeros([2*xs,4*ne]);
% Gi = zeros([2*xs,6*ne]);
% for ll = 1:xs
% for i = 1:2:n-1
% [GWs,HWs] = HGnd(xi(ll),yi(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
%         
%         Gi(2*ll-1:2*ll,3*i-2:3*i+3) = Gi(2*ll-1:2*ll,3*i-2:3*i+3) + GWs;
%         if i == n-1
%             Hi(2*ll-1:2*ll,2*i-1:2*i+2) = Hi(2*ll-1:2*ll,2*i-1:2*i+2) + HWs(1:2,1:4);
%             Hi(2*ll-1:2*ll,1:2) = Hi(2*ll-1:2*ll,1:2) + HWs(1:2,5:6);
%         else
%             Hi(2*ll-1:2*ll,2*i-1:2*i+4) = Hi(2*ll-1:2*ll,2*i-1:2*i+4) + HWs;
%         end
% end
% end
% ui = Gi*f - Hi * u;
% for i = 1:xs
%     uix = ui((i-1)*2+1);
%     uiy = ui(i*2);
%     uia(i) = sqrt(uix^2 + uiy^2);
% end
% 
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
end