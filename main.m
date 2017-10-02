clc
clear all

%% Title



%% Parameters
mu = 0.8e5;
nu = 0.2;
n = 12;
nn = 2*n; %number of unknows
ne = n/2; %number of elements

H = zeros(4*ne);
G = zeros([4*ne,6*ne]);

%% Geometry
x = [-4:2:4,4,4:-2:-4,-4];
y =[-2*ones([1,5]),0,2*ones([1,5]),0];

x(1,end+1) = x(1);
y(1,end+1) = y(1);




for ll = 1:n %collocation point
    for i = 1:2:n-1 % element
        if (ll-i)*(ll-i-1)*(ll-i-2)*(ll-i+n-2) == 0 % local integration
            nodo = ll - i +1;
            if (ll == 1)&&(i == n-1)
                nodo = nodo + n;
            end
            [HWs,GWs] = HGnd(x(ll),y(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
            GWs = Gd(x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu,nodo);
%             HWs = zeros([2,6]);
        else
            [HWs,GWs] = HGnd(x(ll),y(ll),x(i),y(i),x(i+1),y(i+1),x(i+2),y(i+2),nu,mu);
        end
        G(2*ll-1:2*ll,3*i-2:3*i+3) = G(2*ll-1:2*ll,3*i-2:3*i+3) + GWs;
        if i == n-1
            H(2*ll-1:2*ll,2*n-3:2*n) = H(2*ll-1:2*ll,2*n-3:2*n) + HWs(1:2,1:4);
            H(2*ll-1:2*ll,1:2) = H(2*ll-1:2*ll,1:2) + HWs(1:2,5:6);
        else
            H(2*ll-1:2*ll,2*i-1:2*i+4) = H(2*ll-1:2*ll,2*i-1:2*i+4) + HWs;
        end
    end
end

for i = 1:n
    H(2*i-1:2*i,2*i-1:2*i) = zeros(2);
    for j = 1:n
        if i == j
        else
            H(2*i-1:2*i,2*i-1:2*i) = H(2*i-1:2*i,2*i-1:2*i) - H(2*i-1:2*i,2*j-1:2*j);
        end
    end
%     if H(2*i-1,2*i-1)<0
        H(2*i-1:2*i,2*i-1:2*i) = H(2*i-1:2*i,2*i-1:2*i) + [1,0;0,1];
%     end
end


%% Apply BC
A = H;
A(:,[5,17,18]) = [];
A = [A,-(G(:,5)+G(:,7)),-(G(:,23)+G(:,25)),-(G(:,24)+G(:,26))];
f = zeros([36,1]);
f(13) = 1000;
f(17) = -1000;
f(31) = 1000;
f(35) = -1000;
F = G*f;
u = A\F;