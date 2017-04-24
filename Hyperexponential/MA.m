clear all; close all; clc;

% This program computes the steady-state distribution of an M/H_2/N+M
% queue, where we set the service rate mu=1 and the abandonment rate
% alpha>=0.

ctime = cputime;

% number of servers
N = 500;
% N = 5;
% mean patient time
gamma = 1;
% gamma = [5.0, 40.0];
% traffic intensity
rho = 1;
% arrival rate
lambda = N*rho;
% We assume that the transition probabilities A0, A1, and A2 when the level
% is greater than n.
n = N+2000;
% We will compute the steady-state probabilities up to level n+trunc.
trunc = 500;
% parameters to determine the H_2 service time distribution (with mu=1).
p = 0.67407; % c^2=4, r=0.1
m1 = 1/6.7407; % c^2=4, r=0.1
c2 = 4.0;
% p = 0.59155; % c^2=3, r=0.1
% m1 = 1/5.9155; % c^2=3, r=0.1
% c2 = 3.0;



% parameters of the H_2 service time distribution.
q = 1-p;
m2 = (1-p*m1)/q;
mu1 = 1/m1; mu2 = 1/m2;

state = N+1;
epsilon = 1e-6;

K = length(gamma);
pii = zeros(K, n+trunc);
for k = 1:K
    % abandonment rate
    alpha = 1/gamma(k);
    % the R matrices
    R = cell(1, n);
    for i = 1:n
        R{i} = zeros(state, state);
    end
    
    % transition matrices for level n.
    A0 = lambda*eye(state);
    A2 = zeros(state, state);
    A2(1, 1:2) = [N*p*mu1+alpha*(n-1-N), N*q*mu1];
    for m = 2:N
        A2(m, m-1:m+1) = [(m-1)*p*mu2, (N-m+1)*p*mu1+(m-1)*q*mu2+...
            alpha*(n-1-N), (N-m+1)*q*mu1];
    end
    A2(N+1, N:N+1) = [N*p*mu2, N*q*mu2+alpha*(n-1-N)];
    A1 = diag(-lambda-sum(A2, 2));
    
    % compute R{n}.
    G = Gmatrix(A0, A1, A2, epsilon);
    U = A1+A0*G;
    R{n} = -A0/U;
    
    if n == N+1
        A0_1 = [p*lambda*eye(N), zeros(N, 1)]+[zeros(N, 1), q*lambda*eye(N)];
    else
        A0_1 = A0;
    end
    
    % compute R{j} for 1<=j<n.
    for i = 1:n-2
        j = n-i;
        if j < N+2
            A0 = A0_1;
            A0_1 = [p*lambda*eye(j-1), zeros(j-1, 1)]+[zeros(j-1, 1),...
                q*lambda*eye(j-1)];
            A2 = zeros(j, j-1);
            A2(1, 1) = (j-1)*mu1;
            for m = 2:j-1
                A2(m, m-1:m) = [(m-1)*mu2, (j-m)*mu1];
            end
            A2(j, j-1) = (j-1)*mu2;
            if j==2
                A1 = [-(lambda+mu1), 0; 0, -(lambda+mu2)];
            else
                A1 = diag(-lambda-sum(A2, 2));
            end
        else
            A0 = A0_1;
            A0_1 = lambda*eye(state);
            A2 = zeros(state, state);
            A2(1, 1:2) = [N*p*mu1+alpha*(j-1-N), N*q*mu1];
            for m = 2:N
                A2(m, m-1:m+1) = [(m-1)*p*mu2, (N-m+1)*p*mu1+...
                    (m-1)*q*mu2+alpha*(j-1-N), (N-m+1)*q*mu1];
            end
            A2(N+1, N:N+1) = [N*p*mu2, N*q*mu2+alpha*(j-1-N)];
            A1 = diag(-lambda-sum(A2, 2));
        end
        U = A1+A0*G;
        Tmp = -inv(U);
        G = Tmp*A2;
        R{j} = A0_1*Tmp;
    end
    
    % transition probabilities at level 0.
    A0 = [p*lambda, q*lambda];
    A1 = -lambda;
    U = A1+A0*G;
    
    % the steady-state probabilities
    pii(k, :) = piscalar(R, U, n, trunc);
    clear R;
end

ctime = cputime-ctime;
disp(['CPU time used: ', num2str(ctime), ' sec.']);

sig2 = rho+c2+rho-1;
queue = N*(rho-1)*gamma;

save('matrix_n_500_rho_12_H2.mat', 'pii', 'N');

set(0,'DefaultTextInterpreter', 'latex');

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