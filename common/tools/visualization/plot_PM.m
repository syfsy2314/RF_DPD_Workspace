function plot_PM(x,y0,y1)
    %y0是无DPD
    %y1是待比较DPD
    figure;
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
    plot(abs(x), pm_pa_org,'b.' );%黑
    hold on;
    plot(abs(x), pm_pa_mp, 'r.');%红
    hold on;
    title('AM-PM 特性曲线');
    xlabel('归一化输入幅度');
    ylabel('相位变化(rad)');
    legend('No DPD','Propose DPD');
    axis([0,1,-1,1]);
end