clc
clear
load('Result_MVDR.mat')
load('Result_MWF.mat')
SNR_input =[-5,-2,0,2,5];
SNR_mwf = MWF(:,3);
SNR_mvdr = MVDR(:,3);
STOI_mwf = MWF(:,4);
STOI_mvdr = MVDR(:,4);
subplot(2,1,1);
plot(SNR_input,SNR_mvdr,'b-v',SNR_input,SNR_mwf,'r-v');
legend('MVDR','MWF');
xlabel('Relative Input SNR(dB)');
ylabel('Global Output SNR(dB)');
grid on;
subplot(2,1,2);
plot(SNR_input,STOI_mvdr,'b-v',SNR_input,STOI_mwf,'r-v');
legend('MVDR','MWF');
xlabel('Relative Input SNR(dB)');
ylabel('STOI')
grid on;