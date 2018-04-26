% This file is for Asset pricing Assignment 1.
% Yuhao Zhu
% 2014 Jan

clc;

% Import data
start_date = '31122010';
end_date = '31122013'; 
stocks = hist_stock_data(start_date, end_date, 'frequency', 'd', 'tickers.txt'); 
DJIA = hist_stock_data(start_date, end_date, 'DJIA', 'frequency', 'd'); 

% For all components
N=30;
L=length(stocks(1,1).AdjClose);

for i=1:N
    stocks(1,i).R=zeros(L-1,1);
    for j=1:L-1
        stocks(1,i).R(j,1)=(stocks(1,i).AdjClose(j,1)-stocks(1,i).AdjClose(j+1,1))/stocks(1,i).AdjClose(j+1,1);
    end
end

m=zeros(N,1);
for i=1:N
    m(i,1)=mean(stocks(1,i).R);
end

A=zeros(L-1, N);
for i=1:N
    A(:,i)=stocks(1,i).R;
end

V=cov(A);
Std=sqrt(diag(V));

% For DJIA index
DJIA.R=zeros(L-1,1);
for j=1:L-1
    DJIA.R(j,1)=(DJIA(1,1).AdjClose(j,1)-DJIA(1,1).AdjClose(j+1,1))/DJIA(1,1).AdjClose(j+1,1);
end
mDJIA=mean(DJIA.R);
StdDJIA=std(DJIA.R);

% Calculate the EF and plot
e=ones(N,1);
A=e'*(V\m);
B=m'*(V\m);
C=e'*(V\e);
D=B*C-A*A;
mu=-0.0005:0.0001:0.0035;
sigma=sqrt(1/C+C/D*(mu-A/C).^2);
plot(sigma, mu);
axis([0 0.02 -0.0005 0.0035]);
title('Mean-Variance frontier');
ylabel('Mean of portfolio return');
xlabel('Standard deviation of portfolio return');
hold on
scatter(Std, m, 25, 'fill');
scatter(StdDJIA, mDJIA, 75, 'fill');
hold off;

% For any arbitrary
mu_p=0.002;
sigma_p=sqrt(1/C+C/D*(mu_p-A/C).^2);
mu_zp=A/C-D/(C*C)/(mu_p-A/C);
sigma_zp=sqrt(1/C+C/D*(mu_zp-A/C).^2);


% Plot the Zero beta
e=ones(N,1);
A=e'*(V\m);
B=m'*(V\m);
C=e'*(V\e);
D=B*C-A*A;
mu=-0.0005:0.0001:0.0035;
sigma=sqrt(1/C+C/D*(mu-A/C).^2);
plot(sigma, mu);
axis([0 0.02 -0.0005 0.0035]);
title('Mean-Variance frontier and zero-beta portfolio');
ylabel('Mean of portfolio return');
xlabel('Standard deviation of portfolio return');
hold on
scatter(Std, m, 25, 'fill');
scatter(StdDJIA, mDJIA, 75, 'fill');
scatter(sigma_p, mu_p, 75, 'fill');
scatter(sigma_zp, mu_zp, 75, 'fill');
hold off;