clear all; close all; clc;

% This program computes the steady-state distribution of an M/E_2/N+M
% queue, with abandonment rate alpha>0 by directly solving (pi')*G = 0

ctime = cputime;

% number of servers
N = 100;
% N = 5;
% mean patient time
gamma = 1;
% gamma = [5.0, 40.0];
% traffic intensity
rho = 1;
% arrival rate
lambda = N*rho;
%service rate for each phase (have to be identical for E_2)
nu1 = 2;
nu2 = 2;
c2 = 1/2;
% We assume that the transition probabilities A0, A1, and A2 when the level
% is greater than n.
n = N+50;
% We will compute the steady-state probabilities up to level n+trunc.
trunc = 50;
gsize = (n+trunc-N-1)*(N+1) + (N+2)*(N+1)/2; %dimensions of generator matrix
G = sparse(gsize,gsize);
state = N+1;
epsilon = 1e-6;

pii = zeros(1, n+trunc);
% abandonment rate
alpha = 1/gamma;

A0_1 = lambda*eye(state);

A1 = zeros(state, state);
A2 = zeros(state, state);
A2(1, 1:2) = [alpha*(n+trunc-1-N), N*nu2];
for m = 2:N
    A2(m, m:m+1) = [alpha*(n+trunc-1-N), (N-m+1)*nu2];
end
A2(N+1, N+1) = alpha*(n+trunc-1-N);
for m = 2:N+1
    A1(m,m-1) = nu1*(m-1);
end
A1 = A1 - diag(lambda + sum(A2,2) + sum(A1,2));
row_start = (N+2)*(N+1)/2 + (n+trunc-N-2)*(N+1)+1;
row_end = (N+2)*(N+1)/2 + (n+trunc-N-1)*(N+1);
col_start = (N+2)*(N+1)/2 + (n+trunc-N-3)*(N+1)+1;
col_end = (N+2)*(N+1)/2 + (n+trunc-N-1)*(N+1);
G(row_start:row_end, col_start:col_end) = [A2 A1];

% compute R{j} for 1<=j<n.
for i = 1:n+trunc-2
    j = n+trunc-i;
    if j < N+1
        A0 = A0_1;
        A0_1 = [zeros(j-1, 1),lambda*eye(j-1)];
        A1 = zeros(j, j);
        A2 = [nu2*diag(j-1:-1:1); zeros(1,j-1)];
        for m = 2:j
                A1(m,m-1) = nu1*(m-1);
        end
        A1 = A1 - diag(lambda + sum(A2, 2) + sum(A1,2));
        row_start = (j)*(j-1)/2 +1;
        row_end = (j+1)*(j)/2 ;
        col_start = (j-1)*(j-2)/2 +1;
        col_end = (j+2)*(j+1)/2;
    elseif j == N+1
        A0 = A0_1;
        A0_1 = [zeros(j-1, 1),lambda*eye(j-1)];
        A1 = zeros(j, j);
        A2 = [nu2*diag(j-1:-1:1); zeros(1,j-1)];
        for m = 2:j
            A1(m,m-1) = nu1*(m-1);
        end
        A1 = A1 - diag(lambda + sum(A2, 2) + sum(A1,2));
        row_start = (j)*(j-1)/2 +1;
        row_end = (j+1)*(j)/2 ;
        col_start = (j-1)*(j-2)/2 +1;
        col_end = (j+1)*(j)/2 + N+1;
    else
        A0 = A0_1;
        A0_1 = lambda*eye(state);
        A1 = zeros(state, state);
        A2 = zeros(state, state);
        A2(1, 1:2) = [alpha*(j-1-N), N*nu2];
        for m = 2:N
            A2(m, m:m+1) = [alpha*(j-1-N), (N-m+1)*nu2];
        end
        A2(N+1, N+1) = alpha*(j-1-N);
        for m = 2:N+1
            A1(m,m-1) = nu1*(m-1);
        end
        A1 = A1 - diag(lambda + sum(A2, 2) + sum(A1,2));
        row_start = (N+2)*(N+1)/2 + (j-N-2)*(N+1)+1;
        row_end = (N+2)*(N+1)/2 + (j-N-1)*(N+1);
        col_start = (N+2)*(N+1)/2 + (j-N-3)*(N+1)+1;
        col_end = (N+2)*(N+1)/2 + (j-N)*(N+1);
    end
    G(row_start:row_end, col_start:col_end) = [A2 A1 A0];
end

% transition probabilities at level 0.
A0 = [0, lambda];
A1 = -lambda;
G(1,1:3) = [A1 A0];

G(:,end) = 1;
rhs = zeros(1,gsize); 
rhs(end) = 1;
pivec = rhs*inv(G);
ex_pii(1) = pivec(1);
for j = 2:n+trunc
    if j < N+1
        s = j*(j-1)/2+1;
        f = (j+1)*(j)/2;
    else
        s = (N+1)*(N)/2 + (j-N-1)*(N+1) + 1;
        f = (N+1)*(N)/2 + (j-N)*(N+1);
    end 
    ex_pii(j) = sum(pivec(s:f));
end
    

ctime = cputime-ctime;
disp(['CPU time used: ', num2str(ctime), ' sec.']);
    
% sig2 = rho+c2+rho-1;
% queue = N*(rho-1)*gamma;

% save('ma_n_100_beta_0.mat', 'pii', 'N');

% set(0,'DefaultTextInterpreter', 'latex');
% 
% fig1 = figure;
% box on
% xx1 = (-N-queue(1):N*gamma(1)-queue(1))/sqrt(N*gamma(1)*sig2/2);
% line(0:N*gamma(1)+N, pii(1, 1:N*gamma(1)+N+1), 'LineWidth', 1, 'Color', [0.9, 0, 0.1], 'LineStyle', '-');
% hold on
% line(0:N*gamma(1)+N, normpdf(xx1, 0, 1)/sqrt(N*gamma(1)*sig2/2), 'LineWidth', 1, 'Color', [0, 0, 1], 'LineStyle', '-');
% h = legend('matrix-analytic', 'Gaussian');
% set(h, 'Interpreter', 'latex', 'FontSize', 18);
% set(gca, 'FontSize', 12, 'XTick', 50:50:200, 'YTick', 0:0.01:0.03);
% xlim([50, 200]);
% ylim([0, 0.03]);
% xlabel('$X(\infty)$', 'FontSize', 20);
% ylabel('probability', 'FontSize', 20);
% % saveas(fig1, 'MH2100M_gamma1.pdf');
% 
% fig2 = figure;
% box on
% xx2 = (-N-queue(2):N*gamma(2)-queue(2))/sqrt(N*gamma(2)*sig2/2);
% line(0:N*gamma(2)+N, pii(2, 1:N*gamma(2)+N+1), 'LineWidth', 1, 'Color', [0.9, 0, 0.1], 'LineStyle', '-');
% hold on
% line(0:N*gamma(2)+N, normpdf(xx2, 0, 1)/sqrt(N*gamma(2)*sig2/2), 'LineWidth', 1, 'Color', [0, 0, 1], 'LineStyle', '-');
% h = legend('matrix-analytic', 'Gaussian');
% set(h, 'Interpreter', 'latex', 'FontSize', 18);
% set(gca, 'FontSize', 12, 'XTick', 0:200:600, 'YTick', 0:0.002:0.008);
% xlim([0, 600]);
% ylim([0, 0.008]);
% xlabel('$X(\infty)$', 'FontSize', 20);
% ylabel('probability', 'FontSize', 20);
% saveas(fig2, 'MH2100M_gamma10.pdf');

% fig = figure;
% box on
% xx1 = (-N-queue(1):N*gamma(1)-queue(1))/sqrt(N*gamma(1));
% xx2 = (-N-queue(2):N*gamma(2)-queue(2))/sqrt(N*gamma(2));
% line(xx1, pii_den(1, 1:N*gamma(1)+N+1), 'LineWidth', 1, 'Color', [0.9, 0, 0.1], 'LineStyle', '-');
% hold on
% line(xx2, pii_den(2, 1:N*gamma(2)+N+1), 'LineWidth', 1, 'Color', [0.9, 0, 0.1], 'LineStyle', '--');
% hold on
% line(xx2, normpdf(xx2, 0, sqrt(sig2/2)), 'LineWidth', 1, 'Color', [0, 0, 1], 'LineStyle', '-');
% h = legend('matrix-analytic, $\gamma=1.0$',...
%         'matrix-analytic, $\gamma=10$',...
%         'diffusion model');
% set(h, 'Interpreter', 'latex', 'FontSize', 18);
% set(gca, 'FontSize', 12, 'XTick', -6:2:6, 'YTick', 0:0.1:0.3);
% xlim([-6, 6]);
% ylim([0, 0.35]);
% xlabel('the scaled queue length', 'FontSize', 20);
% ylabel('probability density', 'FontSize', 20);
% saveas(fig, 'M_H2_100.pdf');
% fig = figure;
% box on
% xx1 = (-N-queue(1)-1:2*N*gamma(1)-queue(1))/sqrt(N*gamma(1));
% xx2 = (-N-queue(2)-1:2*N*gamma(2)-queue(2))/sqrt(N*gamma(2));
% line(xx1, [0, pii_den(1, 1:2*N*gamma(1)+N+1)], 'LineWidth', 1, 'Color', [0, 0, 1]);
% hold on
% line(xx2, [0, pii_den(2, 1:2*N*gamma(2)+N+1)], 'LineWidth', 1, 'Color', [0, 0.7, 0]);
% hold on
% line(-6:0.05:6, normpdf(-6:0.05:6, 0, sqrt(sig2/2)), 'LineWidth', 1, 'Color', [0.9, 0, 0.1]);
% h = legend('matrix-analytic, $\gamma=5.0$',...
%         'matrix-analytic, $\gamma=40$',...
%         'diffusion model');
% set(h, 'Interpreter', 'latex', 'FontSize', 18);
% set(gca, 'FontSize', 12, 'XTick', -6:2:6, 'YTick', 0:0.1:0.4);
% xlim([-6, 6]);
% ylim([0, 0.4]);
% xlabel('the scaled queue length', 'FontSize', 20);
% ylabel('probability density', 'FontSize', 20);
% saveas(fig, 'M_H2_5.pdf');