clc
clear
load impulse_responses.mat
%% Generate microphone signals
% Target source
N = 577521;
r = 0.5623;
M = 4;
[s1,Fs] = audioread('clean_speech.wav');
if length(s1)<N
    s1=[s1;zeros(N-length(s1),1)];
else
    s1 = s1(1:N);
end
% Interference 1
[n1,~] = audioread('clean_speech_2.wav');
if length(n1)<N
    n1=[n1;zeros(N-length(n1),1)];
else
    n1 = n1(1:N);
end
n1 = n1*r;
% Interference 2
[n2,~] = audioread('babble_noise.wav');
if length(n2)<N
    n2=[n2;zeros(N-length(n2),1)];
else
    n2 = n2(1:N);
end
n2 = n2*r;
% Interference 3
[n3,~] = audioread('aritificial_nonstat_noise.wav');
if length(n3)<N
    n3=[n3;zeros(N-length(n3),1)];
else
    n3 = n3(1:N);
end
n3 = n3*r;
% Interference 4
[n4,~] = audioread('Speech_shaped_noise.wav');
if length(n4)<N
    n4=[n4;zeros(N-length(n4),1)];
else
    n4 = n4(1:N);
end
n4 = n4*r;
% Convoultion to generate 4 microphone signals(noise+speech & pure noise)
X = zeros(4,N+length(h_target)-1);
No = zeros(4,N+length(h_target)-1);
S = zeros(4,N+length(h_target)-1);
for m=1:M
    X(m,:) = conv(h_target(m,:),s1)+conv(h_inter1(m,:),n1)+conv(h_inter2(m,:),n2)...
             +conv(h_inter3(m,:),n3)+conv(h_inter4(m,:),n4);
    No(m,:) = conv(h_inter1(m,:),n1)+conv(h_inter2(m,:),n2)+conv(h_inter3(m,:),n3)...
             +conv(h_inter4(m,:),n4);
    S(m,:) = conv(h_target(m,:),s1);
end
% % sound(X(1,:),Fs);
%% STFT (20ms per frame)
T = 0.02;
N1 = length(X);
K = T*Fs;
L = N1/(Fs*T*0.5)-1;
wind = hann(K).';
% Do STFT for each frame
X_f = zeros(K,L,M);
No_f = zeros(K,L,M);
for m =1:M
    for l=1:L
        X_f(:,l,m) = fft(X(m,0.5*(l-1)*K+1:0.5*(l+1)*K).* wind);
        No_f(:,l,m) = fft(No(m,0.5*(l-1)*K+1:0.5*(l+1)*K).* wind);
    end
end
%% Estimate Rs and ATF & Construct Beamformer per Frequency Bank
W = zeros(M,K);
% Parameter control the signal distortion and noise reduction
mu = 1;
for i=1:K
    % Compute Rx
    xi=squeeze(X_f(i,:,:));
    xi=xi.';
   % xi = xi-mean(xi,2);
    Rxi=xi*xi'/(size(xi,2)-1);
    % Compute Rn
    noi=squeeze(No_f(i,:,:));
    noi=noi.';
   % noi = noi-mean(noi,2);
    Rni=noi*noi'/(size(noi,2)-1);
    % Do GEVD
    [V,D]=eig(Rxi,Rni);
    [eigvals,sidx] = sort(diag(D),'descend');
    U = V(:,sidx);
    Q = inv(U');
    U1 = U(:,1);
    Q1 = Q(:,1);
    sigma = eigvals(1);
    sigma_s = sigma-1;
    e = [1;0;0;0];
    % Optimal filter
    w = (U1*Q1'*e)*sigma_s/(sigma_s+mu);
    W(:,i) = w;
end
%% Reconstruct speech signal and Performance Evaluation ( Based on SNR and STOI with Microphone 1 as the Reference)
S_re_f = zeros(K,L);
S_re = zeros(1,N1);
% Filter each frequency bank per time index
for l=1:L
    for k=1:K
        x = squeeze(X_f(k,l,:));
        S_re_f(k,l)=W(:,k)'* x;
    end
end
% Do IFFT for each time frame and concatenate to form sound signal
for l=1:L
    s_re = ifft(S_re_f(:,l)).';
    S_re(0.5*(l-1)*K+1:0.5*(l+1)*K) = S_re(0.5*(l-1)*K+1:0.5*(l+1)*K)+s_re;
    %S_re((i-1)*k+1:i*k) = s_re;
end

SNR = 10*log10(norm(S(1,:))^2/norm(S_re-S(1,:))^2);
STOI = stoi(S(1,:), S_re, Fs);