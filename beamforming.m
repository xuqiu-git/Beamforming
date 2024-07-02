clear
close all
clc

%% 参数设置
rng default
c=physconst('LightSpeed');      % 光速
fc=1575.42e6;                   % L1频率
fs=5e6;                         % 采样率
lam = c / fc;                   % 波长
load codes_L1CA.mat 
prn = 1;                        % PRN
code0=codes_L1CA(:,1)';         % 伪码
Nc=20;                          % 相干积分时间,ms
N=Nc*fs/1000;                   % 总快拍数
code = code_sample(code0, 1.023e6, fs, N/Nc); % 采样伪码
Nint = 1;                       % 干扰数量
Nsig = 1;                       % 信号数量
incidentAngleI = [85 45 30 15 75;0 0 0 0 0]; % 干扰方向
incidentAngleS = [90;0];        % 信号方向
INR=30;                         % 干噪比
SNR=-25;                        % 信噪比
ant=lam*[-1 -0.5 0 0.5 1; 0 0 0 0 0 ;0 0 0 0 0]; % 阵列位置,5个
kw=2*pi/lam;                    % 波数
de=lam/2;                       % 阵元间距
M=5;                            % 阵元数
array = phased.ConformalArray('ElementPosition',ant);

%% 信号生成
x=zeros(Nsig,N);                % 阵列接收信号
y=ones(Nsig,N);                 % 参考信号
y(1,:)=repmat(code,1,Nc);
yd=despread(y,code);            % 解扩
A=zeros(M,Nsig);                % 导引矢量
gamma=10^(SNR/20);              % 复数幅度
s=gamma*y;                      % 信号生成
Rss = s*s'/N;                   % 信号自相关矩阵
for i=1:Nsig                    % 阵列信号生成
A(:,i)=collectPlaneWave(array,1,incidentAngleS(:,i),fc,c)';
x = x + gamma(i)*A(:,i)*y;
end
%真实相关矩阵
Rxx=x*x'/N;
Ryy=y*y'/N;
Ryx=y*x'/N;
%真实解扩后相关矩阵
xd=despread(x,code);
Rxxd=xd*xd'/Nc;
Ryyd=yd*yd'/Nc;
Ryxd=yd*xd'/Nc;
sd=despread(s,code);
Rssd = sd*sd'/Nc;

%% 参数估计
n = ((randn(M,N))+1j*(randn(M,N)))/sqrt(2); %白噪声
nd=despread(n,code);
xe = x+ n; %接收阵列信号
Q=eye(M); %干扰加噪声协方差矩阵
for i=1:Nint %加干扰
xe = xe+10^(INR/20)*collectPlaneWave(array,1,incidentAngleI(:,i),fc,c)'*(randn(1,N)+1j*randn(1,N))/sqrt(2);
%xe = xe+10^(INR/20)*(collectPlaneWave(array,1,incidentAngleI(:,i),fc,c)'+collectPlaneWave(array,1,[30;0],fc,c)')*(randn(1,N)+1j*randn(1,N))/sqrt(2);
int =  10^(INR/20)*collectPlaneWave(array,1,incidentAngleI(:,i),fc,c)'*(1+1j)/sqrt(2);
Q=Q + int*int';  
end
Rxx=x*x'/N+Q;
Qd=Q*length(code);
xed=despread(xe,code); %解扩阵列接收信号
%ML参数估计
Rxxe=xe*xe'/N;

Ryxe=y*xe'/N;
Ryxed=yd*xed'/Nc;
Rxxed=xed*xed'/Nc;
Qe=Rxxe;

We=sqrtm(inv(Qe));


xed_q=We*xed;
Rxxed_q=xed_q*xed_q'/Nc;

[U,S,V] = svd (Rxxed);
