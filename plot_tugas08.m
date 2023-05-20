clc
clear

% nilai analitik
xa = load('output_xa_pd2_tugas_08.txt');
ya = load('output_ya_pd2_tugas_08.txt');

% nilai numerik
xn = load('output_xn_pd2_tugas_08.txt');
yn = load('output_yn_pd2_tugas_08.txt');

% plot
plot(xa,ya,xn,yn,LineWidth=1.5)
title('Plot Grafik Solusi Numerik dan Solusi Analitik','FontSize',15)
xlabel('nilai x')
ylabel('nilai y')
legend('solusi analitik','solusi numerik','Location','northeast','FontSize',12)
grid on

% Plot the first graph in the first subplot
%subplot(2, 1, 1)
%plot(xa, ya, LineWidth=1.5)
%title('grafik solusi analitik')
%xlabel('x')
%ylabel('y')
%grid on

% Plot the second graph in the second subplot
%subplot(2, 1, 2)
%plot(xn, yn, LineWidth=1.5)
%title('grafik solusi numerik')
%xlabel('x')
%ylabel('y')
%grid on
