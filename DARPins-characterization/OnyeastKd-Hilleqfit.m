EGFP = [1600 800 400 200 100 50 25 12.5 6.25 3.125 1.5625 0.01]; %EGFP concentration range


%3g146 - avg + 2 replicates
Fluorescence1 = [13486 11108 10538 6131 4164 2752 2335 1999 1315 1150 1053 936];  
Fluo1_r1 = [10584 8472 7044 4159 3003 2318 1627 2192 922 1109 978 990]
Fluo1_r2= [16535 13745 14033 8103 5326 3187 3043 1807 1708 1192 1129 883]

%2-g4 - avg + 2 replicates
Fluorescence2 = [13028 11648 9133 8371 9305 6892 5511 3726 2357 1584 1143 856];
Fluo2_r1 = [13393 14926 9417 7628 9863 6686 5577 3984 2514 1493 955 1016]
Fluo2_r2= [12663 8373 8849 9115 8717 7099 5446 3469 2200 1675 1331 696]

%2-h3 - avg + 2 replicates
Fluorescence3 = [33859 39606 32655 25397 27519 20296 22913 18645 13429 10458 6039 623];
Fluo3_r1 = [30027 29824 27687 18944 24282 17931 21245 20383 13973 10900 6557 578]
Fluo3_r2= [37691 49388 37624 31850 30756 23201 24581 16908 12922 10016 5521 668]

lb = [10000 0.5 0];
ub = [45000 10000 1000];
options = optimoptions('lsqcurvefit','ScaleProblem','jacobian');


% MAPPING: Emax = b(1); 42665,  Kd = b(2) background=b(3)
hill_fit = @(b,x)  (b(1).*x./(b(2)+x))+b(3);
b0 = [35000, 200, 0];                                  % Initial Parameter Estimates
[B1,resnorm,residual,exitflag,output] = lsqcurvefit(hill_fit, b0, EGFP, Fluorescence1, lb, ub, options);
AgVct1 = linspace(0.1,10000,1000000);   % Plot Finer Resolution

[B2,resnorm,residual,exitflag,output] = lsqcurvefit(hill_fit, b0, EGFP, Fluorescence2, lb, ub, options);
AgVct2 = linspace(0.1,10000,1000000);   % Plot Finer Resolution

[B3,resnorm,residual,exitflag,output] = lsqcurvefit(hill_fit, b0, EGFP, Fluorescence3, lb, ub, options);
AgVct3 = linspace(0.1,10000,1000000);   % Plot Finer Resolution


Kd = [B1(1,2); B2(1,2); B3(1,2)];

figure('DefaultAxesFontSize',17)
subplot(1,4,1)
plot(EGFP, Fluo1_r1, 'bo')
hold on
plot(EGFP, Fluo1_r2, 'go')
hold on
xlim([1 10000])
ylim([100 100000])
plot(AgVct1, hill_fit(B1,AgVct1), '-r')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off
grid
xlabel('EGFP [nM]')
ylabel('Fluorescence [RFI]')
%legend('Replicate 1', 'Replicate 2', 'Hill Equation Fit', 'Location','SE')
title('3G146')  
ggplotAxes2D([],'AxesTheme','own1','LegendStyle','own1','ColorOrder','Set4');



subplot(1,4,2)
plot(EGFP, Fluo2_r1, 'bo')
hold on
plot(EGFP, Fluo2_r2, 'go')
hold on
xlim([1 10000])
ylim([100 100000])
plot(AgVct2, hill_fit(B2,AgVct2), '-r')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off
grid
xlabel('EGFP [nM]')
ylabel('Fluorescence [RFI]')
%legend('Replicate 1', 'Replicate 2', 'Hill Equation Fit', 'Location','SE')
title('2-G4') 
%legend('Replicate 1', 'Replicate 2', 'Hill Equation Fit', 'Location','SE')
ggplotAxes2D([],'AxesTheme','own1','LegendStyle','own1','ColorOrder','Set4');




subplot(1,4,3)
plot(EGFP, Fluo3_r1, 'bo')
hold on
plot(EGFP, Fluo3_r2, 'go')
hold on
xlim([1 10000])
ylim([100 100000])
plot(AgVct3, hill_fit(B3,AgVct3), '-r')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off
grid
xlabel('EGFP [nM]')
ylabel('Fluorescence [RFI]')
title('2-H3') 
%legend('Replicate 1', 'Replicate 2', 'Hill Equation Fit', 'Location','SE')
ggplotAxes2D([],'AxesTheme','own1','LegendStyle','own1','ColorOrder','Set4');
% 


subplot (1,4,4)
X=categorical({'3G146', '2G-4', '2-H3'})
X=reordercats(X,{'3G146', '2G-4', '2-H3'})
Y=Kd
bar(X,Y)
ylabel('EGFP [nM]')
for i1=1:numel(Y)
    text(X(i1),Y(i1),num2str(Y(i1),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
ggplotAxes2D([],'AxesTheme','own1','LegendStyle','own1','ColorOrder','Set4');
