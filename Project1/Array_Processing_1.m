clc
clear

%% Data generation
d = 2;
M = 5;
N = 20;
Delta = 1/2;
theta = [-20,30]';
f = [0.1,0.12]';
SNR = 20;
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);
s = svd(X);
x = linspace(1,M,M);
% figure(1);
% scatter(x,s','filled');

%% DOA Estimation
theta_0 = esprit(X,d);

%% Frequency Estimation
f_0 = espritfreq(X,d);

%% Joint estimation of directions and frequencies
m=12;
[theta_1,f_1] = joint(X,d,m);

%% Estimation performance
result = zeros(2,4,1000);
Mean = zeros(2,4,6);
Std = zeros(2,4,6);
for i=1:6
    SNR = 4*(i-1);
    for j=1:1000
        [X,A,S] = gendata(M,N,Delta,theta,f,SNR);
        result(:,1,j) = esprit(X,d);
        result(:,2,j) = espritfreq(X,d);
        [result(:,3,j),result(:,4,j)]= joint(X,d,m);
    end
    Mean(:,:,i) = mean(result,3);
    Std(:,:,i) = std(result,0,3);
end
Mean = reshape(Mean,8,6);
Std =  reshape(Std,8,6);

% figure(2);
% x = linspace(0,20,6);
% plot(x,Mean(2,:));
% hold on
% plot(x,Mean(6,:));
% axis([0 20,0 35]);
% legend('esprit','joint');
% xlabel('SNR');
% ylabel('Estimation of angle');
% 
% figure(3);
% plot(x,Mean(3,:));
% hold on
% plot(x,Mean(7,:));
% axis([0 20,-0.2 0.3]);
% legend('esprit','joint');
% xlabel('SNR');
% ylabel('Estimation of frequency');
% 
% figure(4);
% plot(x,Std(2,:));
% hold on
% plot(x,Std(6,:));
% axis([0 20,-5 18]);
% legend('esprit','joint');
% xlabel('SNR');
% ylabel('Standard deviation of angle');
% 
% figure(5);
% plot(x,Std(3,:));
% hold on
% plot(x,Std(7,:));
% axis([0 20,0 0.3]);
% legend('esprit','joint');
% xlabel('SNR');
% ylabel('Standard deviation of frequency');
%% ZF beamformers
SNR=10;
[Xzf,Azf,Szf] = gendata(M,N,Delta,theta,f,SNR);
% With A known
angzfe=esprit(Xzf,d);
Azfe=zeros(M,d);
for i = 1:d
    phizfe = exp(1j*2*pi*Delta*sind(angzfe(i)));
    for k = 0:M-1
       Azfe(k+1,i) = phizfe^k; 
    end
end
W_1 = (pinv(Azfe)).';
S_1 = W_1.'*Xzf;
% With S known
Szfe=zeros(d,N);
freqzfe=espritfreq(Xzf,d);
nzf = (1:N);
for i = 1:d
    Szfe(i,:) = exp(1j*2*pi*nzf*freqzfe(i));
end
W_2 = (pinv(Xzf*Szfe.'*inv(Szfe*Szfe.'))).';
S_2 = W_2.'*Xzf;

% figure(6);
% x = linspace(1,N,N);
% plot(S,'ro');
% hold on
% plot(S_1,'bo');
% hold on
% plot(S_2,'go');
% axis([-1.1 1.1,-1.1 1.1]);
%% compare the spatial response
% angle
respons1=zeros(1,180);
respons2=zeros(1,180);
% first response
for angle=-89:90
the1=exp(1j*pi*sind(angle));
    for i=1:M
        Ath1(i,:)=the1.^(i-1);
    end
respons1(1,angle+90)=dot(W_1'*Ath1,W_1'*Ath1);
end
% second response
for angle=-89:90
the2=exp(1j*pi*sind(angle));
    for i=1:M
        Ath2(i,:)=the2.^(i-1);
    end
respons2(angle+90)=dot(W_2'*Ath2,W_2'*Ath2);
end
figure(7);
plot(-89:90,respons1);
hold on
plot(-89:90,respons2);
legend('ZF1','ZF2');
xlabel('angle');
ylabel('response');
%% Function Definition
function [X,A,S] = gendata(M,N,Delta,theta,f,SNR)
% length(theta) = length(f) = d
n = (1:N);
S = zeros(length(theta),N);
S_n = zeros(length(theta),N);
A = zeros(M,length(theta));
% S & S_n matrix (d x N)
for i = 1:length(theta)
    S(i,:) = exp(1j*2*pi*n*f(i));
    S_n(i,:)= awgn(S(i,:),SNR,'measured');
end
% A matrix (M x d)
for i = 1:length(theta)
    phi = exp(1j*2*pi*Delta*sind(theta(i)));
    for k = 0:M-1
       A(k+1,i) = phi^k; 
    end
end
% X matrix (M x N)
X = A*S_n;
end

function f = esprit(X,d)
[M,~] = size(X);
Z1=X(1:M-1,:);
Z2=X(2:M,:);
Z=[Z1;Z2];
[U,~,~] = svd(Z);
U_z = U(:,1:d);
U_x = U_z(1:M-1,:);
U_y = U_z(M:2*(M-1),:);
f = sort(asind(angle(eig(pinv(U_x)*U_y))/pi));
end

function f = espritfreq(X,d)
[M,N] = size(X);
Z_u = X(:,1:N-1);
Z_l = X(:,2:N);
Z = [Z_u;Z_l];
[U,~,~] = svd(Z);
U_z = U(:,1:d);
U_x = U_z(1:M,:);
U_y = U_z(M+1:2*M,:);
f = sort(angle(eig(pinv(U_x)*U_y))/(2*pi));
end

function [theta,f] = joint(X,d,m)
[M,N] = size(X);
Z = zeros(M*m,N-m+1);
for i=1:N-m+1
    for k=0:m-1
        Z(k*M+1:(k+1)*M,i)=X(:,i+k);
    end
end
[U,~,~] = svd(Z);
U = U(:,1:d);
U_phix = U(1:M*(m-1),:);
U_phiy = U(M+1:M*m,:);
M_phi = pinv(U_phix)*U_phiy;
U_thex = [];
U_they = [];
for i=0:m-1
    U_thex =  vertcat(U_thex,U(i*M+1:(i+1)*M-1,:));
    U_they =  vertcat(U_they,U(i*M+2:(i+1)*M,:));  
end
M_the = pinv(U_thex)*U_they;
% f = angle(eig(M_phi))/(2*pi);
% theta = asind(angle(eig(M_the))/pi);
[T,Theta] = eig(M_the);
theta = sort(asind(angle(eig(Theta))/pi));
Phi = T\M_phi*T;
f = sort(angle(eig(Phi))/(2*pi));
end