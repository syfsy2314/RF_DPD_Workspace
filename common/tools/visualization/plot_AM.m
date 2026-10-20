function plot_AM(x,y0,y1)
    %y0是无DPD
    %y1是待比较DPD
    figure;
    plot(abs(x)/max(abs(x)),abs(y0)/max(abs(y0)),'b.');%黑
    hold on;
    plot(abs(x)/max(abs(x)),abs(y1)/max(abs(y1)),'r.');%红
    hold on;
    title('AM-AM 特性曲线'); 
    xlabel('归一化输入幅度');
    ylabel('归一化输出幅度') ; 
    legend('No DPD','Propose DPD');
    grid on;zoom on; 
end