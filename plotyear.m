function plotyear(DOY,data1,data2,data3)

plot(DOY,data1,'LineWidth',2)
if nargin >= 3
    hold on;
    plot(DOY,data2,'--','LineWidth',2);
    if nargin == 4 
        plot(DOY,data3,'LineWidth',2);
    end
    hold off;
end

    
    xlim([1 365]);grid on;
    set(gca,'XTick',[1,32,60,91,121,152,182,213,244,274,305,335])
    set(gca,'XTickLabel',{'            Jan','            Feb','            Mar',...
        '            Apr','            May','            Jun','            Jul',...
        '            Aug','            Sep','            Oct','            Nov',...
        '            Dec'})
end
