clc
clear all

%% Title



%% Parameters
global mu nu n nn ne
mu = 0.8e5;
nu = 0.2;
n = 28;
nn = 2*n; %number of unknows
ne = n/2; %number of elements

Hma = zeros(4*ne);
Gma = zeros([4*ne,6*ne]);

%% Geometry
x = [0:2:20,20,20:-2:0,zeros([1,5])];
y =[zeros([1,11]),2,4*ones([1,11]),10/3:-2/3:2/3];
xo = x;
yo = y;
x(1,end+1) = x(1);
y(1,end+1) = y(1);


%% BC

f = zeros([3*n,1]);
f(34) = 1500;
f(70) = -7500/9;
f(72) = -12000/9;
f(74) = -12000/9;
f(76) = -1500;
f(78) = -12000/9;
f(80) = -12000/9;
f(82) = -7500/9;
bct = ones([3*n,1]);
dc = [2,67:2:83];
bct(dc) = 0;


%% BUILD H & G


for ll = 1:n %collocation point
    for i = 1:2:n-1 % element
        if (ll-i)*(ll-i-1)*(ll-i-2)*(ll-i+n-2) == 0 % local integration
            nodo = ll - i +1;
            if (ll == 1)&&(i == n-1)
                nodo = nodo + n;
            end
            [~,HWs] = HGnd(x(ll),y(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
            GWs = Gd(x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu,nodo);
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
for i = 1:n
    ux(i) = u(2*i-1);
    uy(i) = u(2*i);
end







node = [xo;yo]';
edge = [1:n;2:n,1]';
hfun = +.5 ; 
[vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun) ;
[nni,~] = size(vert);
ui = zeros(size(vert));
for i = 1:nni
    ui(i,1:2) = InnerDisplacement( x,y,vert(i,1),vert(i,2),ux,uy,n,nu,mu,f,u );
end
figure;
uiabs = sqrt(ui(:,1).^2 + ui(:,2).^2);
    patch('faces',tria(:,1:3),'vertices',vert, ...
'facevertexcdata' , ui(:,2), ...
        'facecolor','interp', ...
        'edgecolor','none') ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;