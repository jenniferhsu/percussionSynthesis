% plot B vs. b0
B = -1:0.01:1;
b0 = (sqrt(1 - B.^2) - 1)./B;
saveDissertationFigures = 1;

figure
plot(B, b0, 'linewidth', 2);
xlabel('B');
ylabel('b_0');
title('Relationship between B and b_0')
grid on
set(gca, 'FontSize', 15)
if saveDissertationFigures==1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, 'figures/BVsb0', 'epsc')
end
