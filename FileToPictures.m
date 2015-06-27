clear; close all;
% load('file_SNR_-15_-15_delNs_0_100_Ntests_1000_Fdop_0Hz.mat');
[FileName,PathNameLux] = uigetfile({'*.mat'},'Choose optical data file.');
%[FileName,PathNameLux] = uigetfile({'*.dat;*.mat'},'Choose optical data file.',pathNameMy);
if FileName ~= 0, %Если есть файл с оптическими данными
    filename = [PathNameLux,FileName]
end
load(filename);
ThreshErr = 1024;%0.5;%0.5;

% FdSchmidlAll_awgn_mean =mean(FdSchmidlAll_awgn,3);
% FdProposedAll_awgn_mean = mean(FdProposedAll_awgn,3);
% 
% FdSchmidlAll_ray_mean =mean(FdSchmidlAll_ray,3);
% FdProposedAll_ray_mean = mean(FdProposedAll_ray,3);

% 
% figure;hold on;grid on; 
% plot(Freqs,FdSchmidlAll_awgn_mean(1,:));
% plot(Freqs,FdProposedAll_awgn_mean(1,:),'r');
% plot(Freqs,FdSchmidlAll_ray_mean(1,:),'g'); 
% plot(Freqs,FdProposedAll_ray_mean(1,:),'k');
% plot(Freqs,Freqs,'m');
% legend('FdSchmidlAll awgn','FdProposedAll awgn','FdSchmidlAll ray','FdProposedAll ray');

FdSchmidlAll_awgnSnr(:,:) = FdMyAwgn(:,1,:);
FdSchmidlAll_raySnr(:,:)=FdMyRay(:,1,:);
FdProposedAll_awgnSnr(:,:) = FdWuAwgn(:,1,:);
FdProposedAll_raySnr(:,:) = FdWuRay(:,1,:);
% FdSchmidlOldAll_awgnSnr(:,:) = FdSchmidlOldAll_awgn(:,1,:);
% FdSchmidlOldAll_raySnr(:,:)=FdSchmidlOldAll_ray(:,1,:);
%ONLY FOR ONE VALUE OF FreqOffset!!!

% [indErrSNawgnRow,indErrSNawgnCol] = find(abs(FdSchmidlAll_awgnSnr-Freqs)<ThreshErr);
% [indErrSNrayRow,indErrSNrayCol] = find(abs(FdSchmidlAll_raySnr-Freqs)<ThreshErr);
% [indErrSOawgnRow,indErrSOawgnCol] = find(abs(FdSchmidlOldAll_awgnSnr-Freqs)<ThreshErr);
% [indErrSOrayRow,indErrSOrayCol] = find(abs(FdSchmidlOldAll_raySnr-Freqs)<ThreshErr);
% [indErrProawgnRow,indErrProawgnCol] = find(abs(FdProposedAll_awgnSnr-Freqs)<ThreshErr);
% [indErrProrayRow,indErrProrayCol] = find(abs(FdProposedAll_raySnr-Freqs)<ThreshErr);

for i=1:length(SNRs)
    ind1 = find(abs(FdSchmidlAll_awgnSnr(i,:)-Freqs)<ThreshErr);
    FdSchmidlAll_awgnSnr1(i) = sqrt(mean(power((FdSchmidlAll_awgnSnr(i,ind1)-Freqs),2)));
    ErRateSchmidlAll_awgnSnr1(i) = (Ntests- length(ind1))/Ntests;
    ind1 = find(abs(FdSchmidlAll_raySnr(i,:)-Freqs)<ThreshErr);
    FdSchmidlAll_raySnr1(i) = sqrt(mean(power((FdSchmidlAll_raySnr(i,ind1)-Freqs),2)));
    ErRateSchmidlAll_raySnr1(i) = (Ntests- length(ind1))/Ntests;
    
    ind1 = find(abs(FdProposedAll_awgnSnr(i,:)-Freqs)<ThreshErr);
    FdProposedAll_awgnSnr1(i) = sqrt(mean(power((FdProposedAll_awgnSnr(i,ind1)-Freqs),2)));
    ErRateProposedAll_awgnSnr1(i) = (Ntests- length(ind1))/Ntests;
    ind1 = find(abs(FdProposedAll_raySnr(i,:)-Freqs)<ThreshErr);
    FdProposedAll_raySnr1(i) = sqrt(mean(power((FdProposedAll_raySnr(i,ind1)-Freqs),2)));
    ErRateProposedAll_raySnr1(i) = (Ntests- length(ind1))/Ntests;
    
%     ind1 = find(abs(FdSchmidlOldAll_awgnSnr(i,:)-Freqs)<ThreshErr);
%     FdSchmidlOldAll_awgnSnr1(i) = sqrt(mean(power((FdSchmidlOldAll_awgnSnr(i,ind1)-Freqs),2)));
%     ErRateSchmidlOldAll_awgnSnr1(i) = (Ntests- length(ind1))/Ntests;
%     ind1 = find(abs(FdSchmidlOldAll_raySnr(i,:)-Freqs)<ThreshErr);
%     FdSchmidlOldAll_raySnr1(i) = sqrt(mean(power((FdSchmidlOldAll_raySnr(i,ind1)-Freqs),2)));
%     ErRateSchmidlOldAll_raySnr1(i) = (Ntests- length(ind1))/Ntests;
end

% FdSchmidlAll_awgnSnr1 = sqrt(mean(power((FdSchmidlAll_awgnSnr-Freqs),2),2));
% FdProposedAll_awgnSnr1 = sqrt(mean(power((FdProposedAll_awgnSnr-Freqs),2),2));
% FdSchmidlAll_raySnr1 = sqrt(mean(power((FdSchmidlAll_raySnr-Freqs),2),2));
% FdProposedAll_raySnr1 = sqrt(mean(power((FdProposedAll_raySnr-Freqs),2),2));
% 
% FdSchmidlOldAll_awgnSnr1 = sqrt(mean(power((FdSchmidlOldAll_awgnSnr-Freqs),2),2));
% FdSchmidlOldAll_raySnr1 = sqrt(mean(power((FdSchmidlOldAll_raySnr-Freqs),2),2));



figure(1);
% plot(SNRs,FdSchmidlAll_awgnSnr1);
% plot(SNRs,FdProposedAll_awgnSnr1,'r');
semilogy(SNRs,FdSchmidlAll_awgnSnr1);hold on;grid on;

semilogy(SNRs,FdProposedAll_awgnSnr1,'r');
% semilogy(SNRs,FdSchmidlOldAll_awgnSnr1,'g');
% legend('Предлагаемый метод','Метод Ву','Исходный Шмидля');
legend('Предлагаемый метод','Метод Ф. Ву');
% legend('Предложенный','Исходный Шмидля');
xlabel('ОСШ, дБ');ylabel('СКО частоты в межканальных интервалах');

figure(2);
% plot(SNRs,FdSchmidlAll_raySnr1);
% plot(SNRs,FdProposedAll_raySnr1,'r');
semilogy(SNRs,FdSchmidlAll_raySnr1);hold on;grid on;

semilogy(SNRs,FdProposedAll_raySnr1,'r');
% semilogy(SNRs,FdSchmidlOldAll_raySnr1,'g');
% legend('Предложенный','Исходный Шмидля');
% legend('Предлагаемый метод','Метод Ву','Исходный Шмидля');
legend('Предлагаемый метод','Метод Ф. Ву');
xlabel('ОСШ, дБ');ylabel('СКО частоты в межканальных интервалах');

% figure(3);
% % plot(SNRs,FdSchmidlAll_awgnSnr1);
% % plot(SNRs,FdProposedAll_awgnSnr1,'r');
% semilogy(SNRs,ErRateSchmidlAll_awgnSnr1);hold on;grid on;
% semilogy(SNRs,ErRateProposedAll_awgnSnr1,'r');
% semilogy(SNRs,ErRateSchmidlOldAll_awgnSnr1,'g');
% legend('Модифицированный Шмидля','Предложенный','Исходный Шмидля');
% xlabel('ОСШ, дБ');ylabel('Оценка вероятности ошибки декодирования');
% 
% 
% figure(4);
% % plot(SNRs,FdSchmidlAll_awgnSnr1);
% % plot(SNRs,FdProposedAll_awgnSnr1,'r');
% semilogy(SNRs,ErRateSchmidlAll_raySnr1);hold on;grid on;
% semilogy(SNRs,ErRateProposedAll_raySnr1,'r');
% semilogy(SNRs,ErRateSchmidlOldAll_raySnr1,'g');
% legend('Модифицированный Шмидля','Предложенный','Исходный Шмидля');
% xlabel('ОСШ, дБ');ylabel('Оценка вероятности ошибки декодирования');








