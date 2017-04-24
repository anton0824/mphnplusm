clear all; close all; % clc;

alpha = 1.0;
beta = 0;%-0.2*sqrt(500);
p = 0.67407;
m1 = 1/6.7407;
cs2 = 4.0;
r = 0.1;
n = 100;

ref = [alpha/(1+cs2) 0.5; -beta/alpha*p, -beta*r;...
    -beta/alpha*(1-p), -beta*(1-r)]; 
degree = 16;
lb1 = -15; ub1 = 40; 
lb2 = -15; ub2 = 40;
% spc = 1.0;
spc = 0.5;
% spc = 0.25;
% spc = 0.125;
numberWithinGrid = spc*40;

mu = 1;
rho = 1-beta/sqrt(n);
rhom = min([rho, 1]);

Gamma = (min(1, rho)+rho)*[p, 0; 0, 1-p];
grid =[lb1:spc:ub1; lb2:spc:ub2];
stp = spc/numberWithinGrid;

starttime = cputime;

% mex lecoeffs.cpp;
[Ar, Ac, Av, b] = lecoeffs(Gamma, grid, degree, alpha, beta, p, m1, ref,n);
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
% mex rpdens.cpp;
outp = rpdens(Gamma, grid, u, numberWithinGrid, alpha, beta, p, m1, ref,n);
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

% ysd = y;
% save('H2ysd_n100_rho1_cs4','ysd');
yc = y;
save('H2yc_n100_rho1_cs4','yc');
load ('matrix_n_100_beta_0_r_01_c_4_alpha_1');

fig2 = figure(2);
box on;
line(xx*sqrt(N)+N, y/sqrt(N), 'LineWidth', 1);
hold on;
line(0:N+500, pii(1, 1:N+501),'LineWidth', 1, 'Color', [1 0 0.2], 'LineStyle', '-');
h = legend('diffusion model', 'matrix-analytic');
set(h, 'Interpreter', 'latex', 'FontSize', 18);
xlabel('number of customers in system', 'FontSize', 20);
ylabel('probability', 'FontSize', 20);
% set(gca, 'FontSize', 12, 'XTick', 50:50:200, 'YTick', 0:0.01:0.03);
% xlim([50, 200]);
% ylim([0, 0.03]);
saveas(fig2, 'MH2100M_gamma1_diffusion.pdf');

% fig = figure(1);
% plot(xx, y, xx, ye, 'LineWidth', 1);
% xlim([-6, 10]);
% ylim([0, 0.3]);

