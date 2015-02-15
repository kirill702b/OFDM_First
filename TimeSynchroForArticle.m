clear all; close all;
N = 1024;
seq1 = 100*ifft(mseq(2,log2(N/2)),N/2)';
%  seq1 = randn(1,N/2)+1i*randn(1,N/2);
seq1 = [randn(1,3*N)+1i*randn(1,3*N),seq1,seq1,randn(1,3*N)+1i*randn(1,3*N)];
% seq1 = awgn([zeros(1,3*N)+1i*zeros(1,3*N),seq1,seq1,zeros(1,3*N)+1i*zeros(1,3*N)],1);
for j=1:length(seq1(1,:))-N
    P1 = sum(conj(seq1(1,j:j+N/2-1)).*seq1(1,j+N/2:j+N/2-1+N/2));
    R1 = sum(abs(conj(seq1(1,j:j+N/2-1)).*seq1(1,j+N/2:j+N/2-1+N/2)));%sum(power(abs(seq1(1,j+N/2:j+N/2-1+N/2)),2));
    RespOfFind(1,j) = power(abs(P1),2)/power(R1,2);
end
%%
close all;
seq2 = 100*ifft(mseq(2,log2(N/4)),N/4)';%
% seq2 = randn(1,N/4)+1i*randn(1,N/4);
seq2 = [randn(1,3*N)+1i*randn(1,3*N),seq2,seq2,-seq2,-seq2,randn(1,3*N)+1i*randn(1,3*N)];
for j=1:length(seq2(1,:))-N
    P1 = sum(conj(seq2(1,j:j+N/4-1)).*seq2(1,j+N/4:j+N/4-1+N/4)) ...
        + sum(conj(seq2(1,j+2*N/4:j+N/4-1+2*N/4)).*seq2(1,j+N/4+2*N/4:j+N/4-1+N/4+2*N/4));
    R1 = sum(abs(conj(seq2(1,j:j+N/4-1)).*seq2(1,j+N/4:j+N/4-1+N/4))) ...
        + sum(abs(conj(seq2(1,j+2*N/4:j+N/4-1+2*N/4)).*seq2(1,j+N/4+2*N/4:j+N/4-1+N/4+2*N/4))) ;% sum(power(abs(seq2(1,j+N/4:j+N/4-1+N/4)),2)) + sum(power(abs(seq2(1,j+N/4+2*N/4:j+N/4-1+N/4+2*N/4)),2));
    RespOfFind(2,j) = power(abs(P1),2)/power(R1,2);
end
% figure;hold on;
% plot(real(seq2));
% plot(imag(seq2),'r');
% 
% figure;plot(RespOfFind(2,:));

%%
seq3 = 100*ifft(mseq(2,log2(N/4)),N/4)';%
% seq3 = randn(1,N/4)+1i*randn(1,N/4);
seq3 = [randn(1,3*N)+1i*randn(1,3*N),seq3,seq3(end:-1:1),conj(seq3),conj(seq3(end:-1:1)),randn(1,3*N)+1i*randn(1,3*N)];
for j=N/2+1:length(seq3(1,:))-N/2
    P1 = 0;R1 = 0;
    for k = 0:N/2
        P1 = P1 + seq3(1,j+k-1)*seq3(1,j-k);
        R1 = R1 + abs(seq3(1,j+k-1)*seq3(1,j-k));%power(abs(seq3(1,j+k-1)),2);
    end  
    RespOfFind(3,j-N/2) = power(abs(P1),2)/power(R1,2);
end
figure;hold on;
plot(real(seq1));
plot(imag(seq1),'r');

figure;plot(RespOfFind(1,:));
figure;plot(RespOfFind(2,:));
figure;plot(RespOfFind(3,:));