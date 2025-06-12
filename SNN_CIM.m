clear,clc,close all
%导入数据
filename = 'G1.txt';
delimiterIn = ' ';
headerlinesIn = 1;
prensent01 = importdata(filename,delimiterIn,headerlinesIn);
dataset = prensent01.data;%读取数据
parameters = prensent01.textdata;
%链接表转化连接矩阵
edges = dataset;
G = graph(edges(:,1),edges(:,2),edges(:,3));
g1 = full(adjacency(G,'weighted'));
g2 = full(adjacency(G));
%链接矩阵转化方程组
tu = importdata(filename);
N = tu(1,1);
bian = -sum(g1,'all');


%参数
d0 = mean(sum(g2,2));
p = 1.0;
C0 = 0;
h = 0.05;
T = 1500;

y = tiledlayout(3,2); %创建2*2画图
%分叉演化曲线
[X,B,R,E,t,C,Z,Eo,b,t1] = Euler(N,g1,C0,d0,h,T,p);%h为步长，T为取值区间
nexttile([1 2])
plot(t,X);%分叉演化曲线
title('振幅-时间演化')

nexttile([1 2])
plot(t,R);%分叉演化曲线
title('纠错-时间演化')

for i = 1:(T/h+1)
Cut(i) = -(bian/2+E(i))/2;
end


nexttile
yyaxis right
plot(t,E);
title('能量-割数演化')
xlabel('t')
ylabel('能量')
hold on
yyaxis left
plot(t,Cut);
ylabel('割数')

% Q = find(B(:,(T/h))>0);
% nexttile%可视化图
% % q=plot(G,'Layout','circle');%环图
% q = plot(G);
% highlight(q,Q,'NodeColor','r')
% title('图结构')

%减小布局四周和每个图块周围的间距
y.Padding = 'compact';
y.TileSpacing = 'compact';
%能量和最大割
format long g
Engry = E(T/h+1);
CutEdge = Cut(T/h+1);
Cuto = -(bian/2+Eo)/2;
Sol = [Engry ,Eo;CutEdge ,Cuto]



%需要的函数定义
function [X,B,R,E,t,C,Z,Eo,b,t1]=Euler(N,g1,C0,d0,h,T,p)%h为步长，T为取值区间，欧拉法求解
ro = 100;
%分配内存
n=round(T/h)+1; %计算离散点的个数
t=zeros(1,n);
X=zeros(N,n);%初始化N个方程为N*n
R=zeros(N,n);
B=zeros(N,n);
E=zeros(1,n);
Z=zeros(N,n);
t1=zeros(1,n);
b=zeros(1,n);
C=C0;
%初始值
X(:,1)=normrnd(0,0.001,[N,1]);%初始化N维方程组的第一列的初始噪声
B(:,1) = X(:,1)./abs(X(:,1));
E(1) = sum(g1*B(:,1).*B(:,1))/2;
R(:,1) = zeros(N,1);
Eo = E(1);
t1(1) = 0.05;
b(:) = 0.3;
for i=2:n            %欧拉法求解振幅
    Z(:,i) = C*(g1*X(:,i-1)); %耦合项
    X(:,i)=X(:,i-1)+h*(p.*X(:,i-1)-X(:,i-1).^3-b(i-1).*R(:,i-1)-tanh(Z(:,i)));%主方程
    R(:,i)=R(:,i-1)+h*(-b(i-1).*R(:,i-1)+X(:,i-1));%主方程
    B(:,i) = X(:,i)./abs(X(:,i));%自旋
    E(:,i) = sum(g1*B(:,i).*B(:,i))/2;%能量
    t1(i) = t1(i-1);
    t(i)=t(i-1)+h;
    dt = t(i)-t1(i);
    C = C+ 0.01*h;
    b(i) = b(i-1)+0.1*h;
    if  E(i) < Eo
           Eo = E(i);
           t1(i) = t(i);
    elseif E(i) > Eo
           b(i) = 0.3; 
    end
    if dt > ro
           t1(i) = t(i);
           b(i) = 0.3;
           C = C0;
    end
end
end