clc
clear all

%% Title



%% Parameters
mu = 0.945e5;
nu = 0.1;
n = 24;
nn = 2*n; %number of unknows
ne = n/2; %number of elements

Hma = zeros(4*ne);
Gma = zeros([4*ne,6*ne]);

%% Geometry
x = 3*cos((-pi/12:-pi/12:-2*pi)-pi/2);
y = 3*sin((-pi/12:-pi/12:-2*pi)-pi/2);




angle = pi/n*2;
nodeangle = (-angle:-angle:-2*pi)-pi/2+angle/2;
nodex = 3*cos(nodeangle);
nodey = 3*sin(nodeangle);
mesh = [1:n;
    2:n,1]';

%% Important Points on Mesh
meshx = nodex(mesh);
meshy = nodey(mesh);
x = ((1*meshx(:,1) + 1*meshx(:,2))/2)';
y = ((1*meshy(:,1) + 1*meshy(:,2))/2)';

x(1,end+1) = x(1);
y(1,end+1) = y(1);
%% BC
f = zeros([3*n,1]);
for i = 1:ne
    f(i*6-5:i*6) = 1.0086*[x(i*2-1)/3*100,y(i*2-1)/3*100,x(i*2)/3*100,y(i*2)/3*100 ...
                    x(i*2+1)/3*100,y(i*2+1)/3*100]';
end
f(33) = 0;
f(52) = 0;
f(69) = 0;
bct = ones([3*n,1]);
bct(33) = 0;
bct(52) = 0;
bct(69) = 0;



for ll = 1:n %collocation point
    for i = 1:2:n-1 % element
        if (ll-i)*(ll-i-1)*(ll-i-2)*(ll-i+n-2) == 0 % local integration
            nodo = ll - i +1;
            if (ll == 1)&&(i == n-1)
                nodo = nodo + n;
            end
            [~,HWs] = HGnd(x(ll),y(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
            GWs = Gd(x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu,nodo);
%             HWs = zeros([2,6]);
        else
            [GWs,HWs] = HGnd(x(ll),y(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
        end
        Gma(2*ll-1:2*ll,3*i-2:3*i+3) = Gma(2*ll-1:2*ll,3*i-2:3*i+3) + GWs;
        if i == n-1
            Hma(2*ll-1:2*ll,2*i-1:2*i+2) = Hma(2*ll-1:2*ll,2*i-1:2*i+2) + HWs(1:2,1:4);
            Hma(2*ll-1:2*ll,1:2) = Hma(2*ll-1:2*ll,1:2) + HWs(1:2,5:6);
        else
            Hma(2*ll-1:2*ll,2*i-1:2*i+4) = Hma(2*ll-1:2*ll,2*i-1:2*i+4) + HWs;
        end
    end
end

for i = 1:n
    Hma(2*i-1:2*i,2*i-1:2*i) = zeros(2);
    for j = 1:n
        if i == j
        else
            Hma(2*i-1:2*i,2*i-1:2*i) = Hma(2cos*i-1:2*i,2*i-1:2*i) - Hma(2*i-1:2*i,2*j-1:2*j);
        end
    end
    if Hma(2*i-1,2*i-1)< 0
        Hma(2*i-1:2*i,2*i-1:2*i) = Hma(2*i-1:2*i,2*i-1:2*i) + [1,0;0,1];
    end
end


%% Apply BC
H = Hma;
G = Gma;
% for i = 1:ne
%     for j = 1:6 % 6 unknowns each element
%         if bct(6*i-6+j) == 0 % 0 is displacement, 1 is traction
%             if (i == ne) && j > 4 % 1st code in the last element
%                 if bct(j-4) > 0
%                     CH = H(:,j-4);
%                     H(:,j-4) =  - G(:,6*i-6+j)*mu;
%                     G(:,6*i-6+j) = - ch;
%                 else
%                     H(:,j-4) = H(:,4*i-4+j)- G(:,6*i-6+j)*mu;
%                     G(:,6*i-6+j) = zeros(size(G(:,6*i-6+j)));
%                 end
%             else
%                 if (i == 1) || (j>2) || (bct(6*i-8+j) ==1)
%                     ch = H(:,4*i-4+j);
%                     H(:,4*i-4+j) =  - G(:,6*i-6+j)*mu;
%                     G(:,6*i-6+j) = - ch;
%                 else
%                     H(:,4*i-4+j) = H(:,4*i-4+j)- G(:,6*i-6+j)*mu;
%                     G(:,6*i-6+j) = zeros(size(G(:,6*i-6+j)));
%                 end
%             end
%         end
%     end
% end
dcG = [52];
dcH = [36];
htemp = Gma(:,dcG);
gtemp = Hma(:,dcH);
Hma(:,dcH) = htemp;
Gma(:,dcG) = gtemp;
                    
%% SOLVE

u = Hma\(Gma*f);