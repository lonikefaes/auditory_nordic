clc;clear all; close all
n = 10;mu = [0 1 2];% artificial data here
X = randn(n,3) + mu;
figure();bp = boxchart(X);hold all;bp.WhiskerLineColor = [1,1,1];% do figure
for it = 1:n
    plot(bp.XData,bp.YData(it,:),'.','MarkerSize',20); % no connecting lines
    plot(bp.XData,bp.YData(it,:),'.-','MarkerSize',20);% lines   
end