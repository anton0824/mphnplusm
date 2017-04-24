clear all; close all; clc;

% This program computes the steady-state distribution of an M/C_2/N+M
% queue, with abandonment rate alpha>0 using the matrix analytic method
%  For documentation about the  M/C_2/N+M model and its
% parameters, see DF_Overloaded.m
% 

ctime = cputime;
tic;
% number of servers
N = 15;
% N = 5;
% mean patient time
gamma = 1;
% gamma = [5.0, 40.0];
% traffic intensity
rho = 1;
%Service Parameters
mean=1;
cs2 = 24; 
nu1 = 2/mean; P12 = 1/(2*cs2); nu2 = nu1*P12;
% arrival rate
lambda = N*rho/mean;
% We assume that the transition probabilities A0, A1, and A2 when the level
% is greater than n.
n = N+2000;
% We will compute the steady-state probabilities up to level n+trunc.
trunc = 500;

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
    A1 = zeros(state, state);
    A2 = zeros(state, state);
    A2(1, 1:2) = [alpha*(n-1-N), N*nu2];
    for m = 2:N
        A2(m, m:m+1) = [alpha*(n-1-N)+(1-P12)*nu1*(m-1), (N-m+1)*nu2];
    end
    A2(N+1, N+1) = alpha*(n-1-N)+(1-P12)*nu1*N;
    for m = 2:N+1
        A1(m,m-1) = nu1*(m-1)*P12;
    end
    A1 = A1 - diag(lambda + sum(A2, 2) + sum(A1,2));
    
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
        if j < N+1
            A0 = A0_1;
            A0_1 = [zeros(j-1, 1),lambda*eye(j-1)];
            A1 = zeros(j, j);
            A2 = zeros(j,j-1);
            A2(1, 1) = nu2*(j-1);
            A2(j,j-1) = (1-P12)*nu1*(j-1);
            for m = 2:j-1
                A2(m, m-1:m) = [(1-P12)*nu1*(m-1), nu2*(j-1-m+1)];
            end
            for m = 2:j
                    A1(m,m-1) = nu1*(m-1)*P12;
            end
            A1 = A1 - diag(lambda + sum(A2, 2) + sum(A1,2));
        elseif j == N+1
            A0 = A0_1;
            A0_1 = [zeros(j-1, 1),lambda*eye(j-1)];
            A1 = zeros(j, j);
            A2 = zeros(j,j-1);
            A2(1, 1) = nu2*(j-1);
            A2(j,j-1) = (1-P12)*nu1*(j-1);
            for m = 2:j-1
                A2(m, m-1:m) = [(1-P12)*nu1*(m-1), nu2*(j-1-m+1)];
            end
            for m = 2:j
                A1(m,m-1) = nu1*(m-1)*P12;
            end
            A1 = A1 - diag(lambda + sum(A2, 2) + sum(A1,2));
        else
            A0 = A0_1;
            A0_1 = lambda*eye(state);
            A1 = zeros(state, state);
            A2 = zeros(state, state);
            A2(1, 1:2) = [alpha*(j-1-N), N*nu2];
            for m = 2:N
                A2(m, m:m+1) = [alpha*(j-1-N)+(1-P12)*nu1*(m-1), (N-m+1)*nu2];
            end
            A2(N+1, N+1) = alpha*(j-1-N)+(1-P12)*nu1*N;
            for m = 2:N+1
                A1(m,m-1) = nu1*(m-1)*P12;
            end
            A1 = A1 - diag(lambda + sum(A2, 2) + sum(A1,2));
        end
        U = A1+A0*G;
        Tmp = -inv(U);
        G = Tmp*A2;
        R{j} = A0_1*Tmp;
    end
    
    % transition probabilities at level 0.
    A0 = [0, lambda];
    A1 = -lambda;
    U = A1+A0*G;
    
    % the steady-state probabilities
    pii(k, :) = piscalar(R, U, n, trunc);
    clear R;
end

toc;
ctime = cputime-ctime;
disp(['CPU time used: ', num2str(ctime), ' sec.']);


save('maC2_n_15_rho_1_cs_24.mat', 'pii', 'N'); 
