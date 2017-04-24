clear all; close all; % clc;

%This code computes the diffusion approximation for the stationary 
%distribution of the customer count in 
%the M/C_2/n+M system. This system is an n-server queue with a 2-phase 
%Coxian service time distribution. 
%The algorithm used is based on the
%procedure in Dai-He, Stochastic Systems, 2013.

%The parameters that can be modified are 
%rho -- the utilization of the system
%cs2 -- squared coefficient of variation of the service time distribution
%alpha -- abandonment rate of the system
%n -- number of servers

%The output is a vector y, which contains the density of the diffusion
%approximation of the total customer count. 

cs2 = 24; %can modify
rho = 1; %can modify
alpha = 1.0; %can modify 
n = 15; %number of servers

%A 2-phase Coxian distribution has parameters nu1, nu2, P12.
%Phases 1 and 2 serve at rates nu1 and nu2. After phase 1, move on to
%phase 2 with probability P12, or leave system with probability 1-P12.
%A common choice of nu1,nu2,P12 is to first pick the desired mean m and 
%coefficient of squared variation cs2, and then set
%nu1 = 2/m, P12=1/(2*cs2), nu2 = P12*nu1
m=1; %mean of the service time distribution. Always set to one.
nu1 = 2/m; P12 = 1/(2*cs2); nu2 = nu1*P12;
r = (1/nu1)/m; %This is utilization rate for server 1
p=1; %always start in phase 1
beta = (1-rho)*sqrt(n); %Halfin-Whitt parameter

%Parameters for the reference density
ref = [alpha/(1+cs2) 0.5; -beta/alpha*p, -beta*r;...
    -beta/alpha*(1-p), -beta*(1-r)]; 
degree = 16;
%Integration limits
lb1 = -10; ub1 = 35; 
lb2 = -10; ub2 = 35;

%Mesh resolution
% spc = 1.0;
spc = 0.5;
% spc = 0.25;
% spc = 0.125;
numberWithinGrid = spc*40;

mu = 1/m;
lambda = n*rho*mu;
del = 1/sqrt(n);

Gamma = [0 0;0 0];
grid =[lb1:spc:ub1; lb2:spc:ub2];
stp = spc/numberWithinGrid;

starttime = cputime;

% The function lecoeffs sets up the system of equations in (3.14) of
% Dai-He 2013. The function is implemented in c++, and the code is
%  contained in lecoeffs.cpp. 

%Must run the following mex command anytime the cpp file is modified.
% mex lecoeffs.cpp;
[Ar, Ac, Av, b] = lecoeffs(Gamma, grid, degree, alpha, lambda, nu1, nu2, P12, n, ref);
M = length(Ar)/36;
nzl = length(nonzeros(Ar));
A = sparse(Ar(1:nzl), Ac(1:nzl), Av(1:nzl), M, M);

disp(['CPU time used for constructing A and b: ', num2str(cputime-starttime),...
    ' sec.']);

invstarttime = cputime;
% disp(['The condition number of A is ', num2str(condest(A), '%10.2e')]);
u = A\b;

% u_p = zeros(length(b), 1);
% [c,R,P] = qr(A, b);
% rr = sum(diag(R)~=0);
% u_p(1:rr) = R(1:rr, 1:rr)\c(1:rr);
% u = P*u_p;

% u_p = zeros(length(b), 1);
% u = u_p;
% [c, R, pm] = qr(A, b, 0);
% rr = sum(diag(R)~=0);
% u_p(1:rr) = R(1:rr, 1:rr)\c(1:rr);
% u(pm) = u_p;


disp(['CPU time used for inverting A: ', num2str(cputime-invstarttime),...
    ' sec.']);
disp(['CPU time used: ', num2str(cputime-starttime), ' sec.']);

clear Ar Ac Av A b;

densstarttime = cputime;

%Compute the diffusion density based on the display below (3.15) in p. 113 
%of Dai-He 2013


% mex rpdens.cpp;
outp = rpdens(Gamma, grid, u, numberWithinGrid, alpha, lambda, nu1, nu2, P12, n, ref);
dimesh = sqrt(length(outp(:, 3)));
Prob = reshape(outp(:, 3), dimesh, dimesh);
clear outp;
Prob = Prob/sum(sum(Prob.*stp^2));
disp(['CPU time used for constructing density: ',...
    num2str(cputime-densstarttime), ' sec.']);
disp(['Total CPU time used for constructing density: ',...
    num2str(cputime-starttime), ' sec.']);

dim = (size(grid, 2)-1)*numberWithinGrid+1;
Prl = fliplr(Prob);
clear Prob

set(0,'DefaultTextInterpreter', 'latex');

y = zeros(1, 2*dim-1);
xx = (lb1+lb2):stp:(ub1+ub2);

% compute the one-dimensional stationary density
for m = 1:2*dim-1
    y(m) = sum(diag(Prl, dim-m)*stp^2)/stp;
end


ysd = y;
save('C2ysd_n15_rho1_cs24','ysd');
% yc = y;
% save('C2yc_n15_rho1_cs24','yc');
load ('maC2_n_15_rho_1_cs_24');

fig2 = figure(2);
box on;
line(xx*sqrt(N)+N, y/sqrt(N), 'LineWidth', 1, 'Color', [1 0 0]);
hold on;
line(0:N+500, pii(1, 1:N+501),'LineWidth', 1, 'Color', [0 0 1], 'LineStyle', '-');
h = legend('state dep.', 'const.', 'matrix-analytic');
set(h, 'Interpreter', 'latex', 'FontSize', 18);
xlabel('number of customers in system', 'FontSize', 20);
ylabel('probability', 'FontSize', 20);

