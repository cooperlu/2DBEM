clc
clear all

%% Title



%% Parameters
mu = 0.945e5;
nu = 0.1;
n = 48;
nn = 2*n; %number of unknows
ne = n/2; %number of elements

Hma = zeros(4*ne);
Gma = zeros([4*ne,6*ne]);

%% Geometry
x = [0:0.2:4,4,4,4,4:-0.2:0,0,0,0];
y =20.*[zeros([1,21]),0.05:0.05:0.15,0.2*ones([1,21]),0.15:-0.05:0.05];

x(1,end+1) = x(1);
y(1,end+1) = y(1);


%% BC

f = zeros([3*n,1]);
f(61:2:71) = 1000;

bct = ones([3*n,1]);
dc = [133:144];
bct(dc) = 0;

%% BUILD H & G


for ll = 1:n %collocation point
    for i = 1:2:n-1 % element
        if (ll-i)*(ll-i-1)*(ll-i-2)*(ll-i+n-2) == 0 % local integration
            nodo = ll - i +1;
            if (ll == 1)&&(i == n-1)
                nodo = nodo + n;
            end
            [~,HWs] = HGndm(x(ll),y(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
            GWs = Gd(x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu,nodo);
        else
            [GWs,HWs] = HGndm(x(ll),y(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
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
            Hma(2*i-1:2*i,2*i-1:2*i) = Hma(2*i-1:2*i,2*i-1:2*i) - Hma(2*i-1:2*i,2*j-1:2*j);
        end
    end
    if Hma(2*i-1,2*i-1)< 0
        Hma(2*i-1:2*i,2*i-1:2*i) = Hma(2*i-1:2*i,2*i-1:2*i) + [1,0;0,1];
    end
end



%% Apply BC
H = Hma;
G = Gma;

for i = 1:ne
    for j = 1:6 % 6 unknowns each element
        if bct(6*i-6+j) == 0 % 0 is displacement, 1 is traction
            if (i == ne) && j > 4 % 1st code in the last element
                if bct(j-4) > 0
                    ch = H(:,j-4);
                    H(:,j-4) =  - G(:,6*i-6+j)*mu;
                    G(:,6*i-6+j) = - ch;
                else
                    H(:,j-4) = H(:,4*i-4+j)- G(:,6*i-6+j)*mu;
                    G(:,6*i-6+j) = zeros(size(G(:,6*i-6+j)));

                end
            else
                if (i == 1) || (j>2) || (bct(6*i-8+j) ==1)
                    ch = H(:,4*i-4+j);
                    H(:,4*i-4+j) =  - G(:,6*i-6+j)*mu;
                    G(:,6*i-6+j) = - ch;
                else
                    H(:,4*i-4+j) = H(:,4*i-4+j)- G(:,6*i-6+j)*mu;
                    G(:,6*i-6+j) = zeros(size(G(:,6*i-6+j)));
                end
            end
        end
    end
end
                  
%% SOLVE

u = H\(G*f*100);

%% Organize


u = u./100;
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


for i = 1:19
    for j = 1:79
        xi((i-1)*79+j) = j*0.05;
        yi((i-1)*79+j) = i*0.01;
    end
end


[~,xs] = size(xi);

Hi = zeros([2*xs,4*ne]);
Gi = zeros([2*xs,6*ne]);
for ll = 1:xs
for i = 1:2:n-1
[GWs,HWs] = HGnd(xi(ll),yi(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
        
        Gi(2*ll-1:2*ll,3*i-2:3*i+3) = Gi(2*ll-1:2*ll,3*i-2:3*i+3) + GWs;
        if i == n-1
            Hi(2*ll-1:2*ll,2*i-1:2*i+2) = Hi(2*ll-1:2*ll,2*i-1:2*i+2) + HWs(1:2,1:4);
            Hi(2*ll-1:2*ll,1:2) = Hi(2*ll-1:2*ll,1:2) + HWs(1:2,5:6);
        else
            Hi(2*ll-1:2*ll,2*i-1:2*i+4) = Hi(2*ll-1:2*ll,2*i-1:2*i+4) + HWs;
        end
end
end
ui = Gi*f - Hi * u;
for i = 1:xs
    uix = ui((i-1)*2+1);
    uiy = ui(i*2);
    uia(i) = sqrt(uix^2 + uiy^2);
end

for i = 1:19
    for j = 1:79
        xim(i,j) = xi((i-1)*79+j);
        yim(i,j) = yi((i-1)*79+j);
        uim(i,j) = uia((i-1)*79+j);
    end
end

% surf(xim,yim,uim)
ux = 0;
uy = 0;
for i = 1:n
    ux(i) = u(2*i-1);
    uy(i) = u(2*i);
end
ux(end+1) = ux(1);
uy(end+1) = uy(1);
coe = 1e2;
plot(x,y,x+ux*coe,y+uy*coe);
legend('original','1000x deformation');