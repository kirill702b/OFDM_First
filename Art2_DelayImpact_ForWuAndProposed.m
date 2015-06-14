clear; close all;
SNR = 50;
dev = 50;
Nfft=2^17;
FreqOffset = 20.4;
N=1024;%number of bits in preamble
EnableGraphs = 0;
CP = 0.1; % Cyclic prefix length as part of symbol period
Nnoise = 20;
FreqDop = 0;
delN = 0;
frs3 = (randi(2,[1 N])-1.5)*2;% вторая последовательность
BP = ifft(frs3,N);

tic;
seq3UpFcAwgnF0Dec500LpfFs1=SequenceThroughChannelFunction(BP,N,FreqOffset,SNR,0,Nnoise,FreqDop,delN);
toc;
if EnableGraphs
    figure;plot(abs(seq3UpFcAwgnF0Dec500LpfFs1));hold on;plot(Nnoise+(1:length(BP)),abs(BP),'r');
end
%% Proposed
maxAwgnIndPro = Nnoise+1;
seq12df(1,:) = seq3UpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndPro:maxAwgnIndPro-1+N);%АГБШ
FdWu = WuMethod(seq12df(1,1:N),BP)
FdMy = MyMethod(seq12df(1,1:N),BP,Nfft,dev)
