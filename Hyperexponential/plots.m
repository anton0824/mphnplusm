clear all; close all; % clc;

lb1 = -15; ub1 = 40;
lb2 = -15; ub2 = 40;
stp=1/40;
xx = (lb1+lb2):stp:(ub1+ub2);


nvals = 100%[10 20 40 80];
cerrs = zeros(1, length(nvals));
sderrs = zeros(1, length(nvals));
kcerrs = zeros(1, length(nvals));
ksderrs = zeros(1, length(nvals));
exacts = zeros(1, length(nvals));
exact_tails = zeros(1, length(nvals));
fig3 = figure(3);
for count = 1:length(nvals)
    n = nvals(count);
    f1 = strcat('H2yc_n', int2str(n) , '_rho1_cs4');
    f2 = strcat('H2ysd_n', int2str(n) , '_rho1_cs4');
    f3 = strcat('matrix_n_', int2str(n) , '_beta_0_r_01_c_4_alpha_1');
    load (f1);
    load (f2);
    load (f3);

    yc(yc<0)=0; ysd(ysd<0)=0;

    fig2 = figure(2);
    box on;
    line(0:N+500, pii(1, 1:N+501),'LineWidth', 1, 'Color', [0 0 1], 'LineStyle', '-');
    line(xx*sqrt(N)+N, yc/sqrt(N), 'LineWidth', 1, 'Color', [0 0.3 0]);
    line(xx*sqrt(N)+N, ysd/sqrt(N), 'LineWidth', 1, 'Color', [1 0 0]);
    h = legend('exact', 'constant','state-dependent');
    hold on;
    set(h, 'Interpreter', 'latex', 'FontSize', 18);
    xlabel('number of customers in system', 'FontSize', 20);
    ylabel('probability', 'FontSize', 20);
    title('Steady-state total customer count distribution', 'FontSize', 26);
    % set(gca, 'FontSize', 12, 'XTick', 50:50:200, 'YTick', 0:0.01:0.03);
    % xlim([50, 200]);
    % ylim([0, 0.03]);
    % saveas(fig2, 'MH2100M_gamma1_diffusion.pdf');

    % fig = figure(1);
    % plot(xx, y, xx, ye, 'LineWidth', 1);
    % xlim([-6, 10]);
    % ylim([0, 0.3]);
    TOL = 0.001;
    clear rel_const
    clear rel_statedep
    pii = pii(1,:);
    exact_TCDF = 1-cumsum(pii);
    const_TCDF = 1-cumsum(yc)*stp;
    statedep_TCDF = 1-cumsum(ysd)*stp;
    UB=sum(exact_TCDF > TOL);
    abs_const = zeros(1,UB); abs_statedep = zeros(1,UB);
    rel_const = zeros(1,UB); rel_statedep = zeros(1,UB);
    as = xx*sqrt(N)+N;
    for l=1:UB
        ind = sum(as < l)+1;
        tail_exact = sum(pii(l+1:end));
        tail_const = sum(yc(ind:end))*stp + (as(ind)-l)/sqrt(N)*yc(ind);
        tail_statedep = sum(ysd(ind:end))*stp + (as(ind)-l)/sqrt(N)*ysd(ind);
        abs_const(l)=abs(tail_const -tail_exact);
        abs_statedep(l) = abs(tail_statedep-tail_exact);
        rel_const(l)=abs(tail_const/tail_exact-1);
        rel_statedep(l) = abs(tail_statedep/tail_exact-1);
    end
    exact_tails(count) = tail_exact;
    kcerrs(count) = rel_const(UB);
    ksderrs(count) = rel_statedep(UB);
    figure(3)
%     subplot(2,4,count);
    plot(exact_TCDF(1:UB),rel_const,'.',exact_TCDF(1:UB),rel_statedep, 'o')
    set(gca,'XDir','reverse')
    xlabel('\fontsize{18}{0}\selectfont$P(T \geq k)$ ','Interpreter','LaTex')
    
    
    figure(4)
    plot(exact_TCDF(1:UB),abs_const,'.',exact_TCDF(1:UB),abs_statedep, 'o')
    set(gca,'XDir','reverse')
    xlabel('\fontsize{18}{0}\selectfont$P(T \geq k)$ ','Interpreter','LaTex')
    
%     ylabh = get(gca,'YLabel');
%     ylabel('$\Big|\frac{P(W \geq x)}{P(T \geq x)} - 1\Big| $','Interpreter','LaTex', 'rot',0)
%     set(ylabh,'Position',get(ylabh,'Position') + [.02 0 0])
    title(strcat('n=',int2str(n)), 'FontSize', 26)
    %Second moment

    %compute exact 2nd moment
    mom = 1;
    ls = ((0:length(pii)-1) - N)/sqrt(N);
    eX2 = sum((abs(ls).^mom).*pii);

    %constant coefficient
    cX2 = sum((abs(xx).^mom).*yc)*stp;

    %state dependent 
    sdX2 = sum((abs(xx).^mom).*ysd)*stp;
    
    exacts(count) = eX2;
    cerrs(count) = abs(cX2-eX2);
    sderrs(count) = abs(eX2-sdX2);
%     [eX2 cX2 sdX2; 0 abs(cX2-eX2) abs(eX2-sdX2)]
end
