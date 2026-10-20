function plot_AM_PM(x,y0,y_dpd)
%x是pa输入
%y0是无dpd的pa输出
%y_dpd是dpd之后的信号
    y1 = y_dpd;
    %subplot(1,2,1);
    yyaxis right
    plot(abs(x)/max(abs(x)),abs(y0)/max(abs(y0)),'k.');%黑
    hold on;
    plot(abs(x)/max(abs(x)),abs(y1)/max(abs(y1)),'r.');%红
    hold on;
    title('');
    xlabel('Normalized Input Magnitude');
    ylabel('Normalized Output Magnitude') ;
    %legend('without DPD','ILC-FDE','Location','northwest');
    set(gca, 'ycolor', [0 0 0]);
    grid on;zoom on;
    %subplot(1,2,2);
    yyaxis left
    PeakReal = max(abs(real(x)));
    PeakImag = max(abs(imag(x)));
    Peak = max([PeakReal PeakImag]);
    x = double(x./Peak);
    y0 = double(y0./Peak);
    y1 = double(y1./Peak);
    pm_pa_org = angle(y0) - angle(x);
    pm_pa_org(pm_pa_org < -pi) = pm_pa_org(pm_pa_org < -pi) + 2 * pi;
    pm_pa_org(pm_pa_org > pi) = pm_pa_org(pm_pa_org > pi) - 2 * pi;
    pm_pa_mp = angle(y1) - angle(x);
    pm_pa_mp(pm_pa_mp < -pi) = pm_pa_mp(pm_pa_mp < -pi) + 2 * pi;
    pm_pa_mp(pm_pa_mp > pi) = pm_pa_mp(pm_pa_mp > pi) - 2 * pi;
    grid on;hold on;
    plot(abs(x), pm_pa_org,'k.' );%黑
    hold on;
    plot(abs(x), pm_pa_mp, 'r.');%红
    hold on;
    title('');
    xlabel('Normalized Input Magnitude');
    ylabel('Phase Shift(rad)');
    set(gca, 'ycolor', [0 0 0]);
    legend('W/o DPD','LD-PLS-DPD','Location','southeast');
    axis([0,1,-1,1]);
end

