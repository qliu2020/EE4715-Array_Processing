clc
clear
%% Signal model
N = 500;
P = 4;
sigma = 0.5;
ind = [-3,-1,1,3];
k = randi(4,N,1);
s = zeros(1,N);
for i=1:N
    s(i) = exp(1j*pi/4*ind(k(i)));
end
x = gendata_conv(s,P,N,sigma);

%% Construct X
X = zeros(2*P,N-1);
for i=1:N-1
    X(:,i)=x((i-1)*P+1:(i+1)*P);
end
R = rank(X);
fprintf('The rank of X is %d.\n',R);

%% Double P
% P_1 = P*2;
% x_1 = gendata_conv(s,P_1,N,sigma);
% X_1 = zeros(2*P_1,N-1);
% for i=1:N-1
%     X_1(:,i)=x_1((i-1)*P_1+1:(i+1)*P_1);
% end
% disp(rank(X_1));

%% Zero-forcing and Wiener equalizer
h = zeros(P,1);
for i=1:P
    m = (i-1)/P;
    if (m>=0 && m<0.25) || (m>=0.5 && m<0.75)
            h(i) = 1;        
    elseif (m>=0.25 && m<0.5) || (m>=0.75 && m<1)
            h(i) = -1;
    else
            h(i) = 0;
    end           
end
H = [zeros(P,1),h;h,zeros(P,1)];
S = [s(2:N);s(1:N-1)];

W_zf = pinv(H)';
id = eye(2*P,2*P);
W_wn = (H*H'+(sigma^2)*id)\H;
disp(W_zf);
disp(W_wn)
lambda_zf=[];
for i=1:2
    f = (H'* W_zf(:,i))/norm(W_zf(:,i));
    fi = abs(f(i)) - abs(f(2/i));
    lambda_zf = cat(2,lambda_zf,fi);
end
disp(lambda_zf);   

lambda_wn=[];
for i=1:2
    f = (H'* W_wn(:,i))/norm(W_wn(:,i));
    fi = abs(f(i)) - abs(f(2/i));
    lambda_wn = cat(2,lambda_wn,fi);
end
disp(lambda_wn);
S_zf = W_zf'*X;
S_wn = W_wn'*X;

real_zf = real(S_zf(1,:));
imag_zf= imag(S_zf(1,:));

real_wn = real(S_wn(1,:));
imag_wn = imag(S_wn(1,:));

real_gt = real(s);
imag_gt = imag(s);

plot(real_zf,imag_zf,'g.',real_wn,imag_wn,'b.',real_gt,imag_gt,'r*');
hold on
xlim = get(gca,'Xlim');
ylim = get(gca,'Ylim');
plot(xlim,[0,0],'k--');
plot([0,0],ylim,'k--');
hold off
legend('ZF','WN','GT')
%% Function defination
function x = gendata_conv(s,P,N,sigma)
x = zeros(N*P,1);
n_dB = 10*log10(sigma^2);
n = wgn(N*P,1,n_dB,'complex');
for t=1:N*P
    for k=1:N
        i = (t-1)/P-(k-1);
        if (i>=0 && i<0.25) || (i>=0.5 && i<0.75)
                h = 1;
        elseif (i>=0.25 && i<0.5) || (i>=0.75 && i<1)
                h = -1;
        else
                h = 0;
        end
        x(t) = x(t)+h*s(k); 
    end
    x(t) = x(t)+n(t);
end
end