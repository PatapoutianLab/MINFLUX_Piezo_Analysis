
% This code analyzes the distance ratios between each blade in a PIEZO1 trimer
% Written by Rachel Lee and Eric Mulhall

clc
clear
close all

%%Test file 1
files2test = {['.mat']}; 

%%Test file 2
files2test = {['.mat']}; 

%%Test file 3
files2test = {['.mat']}; 

%%Save name 1
saveNames = {''};
%%Save name 2
saveNames = {''};
%%Save name 3
saveNames = {''};

cmap = [0 0 0 ; ... Black
        230 159 0 ; ... Orange
        86 180 233 ; ... Sky blue
        0 158 115 ; ... Blueish green
        240 228 66 ; ... Yellow
        0 114 178 ; ... Blue
        213 94 0 ; ... Vermillion
        204 121 167]/255; % Reddish purple

for kk = 1%:length(files2test)

    load(files2test{kk})
    if kk == 33
        kk = 3;
    end

    L = NaN*ones(length(particles),1);
    loc3 = cell(length(particles),1);
    distEach = NaN*ones(length(particles),3);
    distSort = distEach;
    for jj = 1:length(particles)
        L(jj) = size(particles{jj}.coords,1);

        if L(jj) == 3
            loc3{jj} = particles{jj}.coords(:,1:3);

            distEach(jj,:) = pdist(loc3{jj});
            distSort(jj,:) = sort(distEach(jj,:));
        end
    end

    %% Distance Distributions

   fig = figure(1)
   h = histogram(distSort(:,1),10:5:60,'FaceColor',cmap(3,:),'LineWidth',1)

    set(gca,'FontSize',24,'LineWidth',1,'FontName','Arial')
    xlabel('Distance (nm)','FontSize',24,'FontName','Arial')
    ylabel('Count','FontSize',24,'FontName','Arial')
    title(["Shortest", "Inter-Blade Distance"])

    % auto-set axes based on maximum bin height + 25%
    bin_counts1 = get(h, 'Values')';
    max_count = max(bin_counts1);
    ylim([0,max_count*1.25])

    box off
    saveas(gcf,[saveNames{kk} '_ShortEdge.pdf'],'pdf')
%%

   fig = figure(2)
   h = histogram(distSort(:,2),10:5:60,'FaceColor',cmap(4,:),'LineWidth',1)

    set(gca,'FontSize',24,'LineWidth',1,'FontName','Arial')
    xlabel('Distance (nm)','FontSize',24,'FontName','Arial')
    ylabel('Count','FontSize',24,'FontName','Arial')
    title(["Middle", "Inter-Blade Distance"])

    % auto-set axes based on maximum bin height + 25%
    bin_counts2 = get(h, 'Values')';
    max_count = max(bin_counts2);
    ylim([0,max_count*1.25])

    box off
    saveas(gcf,[saveNames{kk} '_MiddleEdge.pdf'],'pdf')
%%

   fig = figure(3)
   h =  histogram(distSort(:,3),10:5:60,'FaceColor',cmap(2,:),'LineWidth',1)

    set(gca,'FontSize',24,'LineWidth',1,'FontName','Arial')
    xlabel('Distance (nm)','FontSize',24,'FontName','Arial')
    ylabel('Count','FontSize',24,'FontName','Arial')
    title(["Longest", "Inter-Blade Distance"])

    % auto-set axes based on maximum bin height + 25%
    bin_counts3 = get(h, 'Values')';
    max_count = max(bin_counts3);
    ylim([0,max_count*1.25])

    box off
    saveas(gcf,[saveNames{kk} '_LongEdge.pdf'],'pdf')

    %% Ratio Distributions

    figure(4)
    h = histogram(distSort(:,1)./distSort(:,2),0:0.1:1.2,'FaceColor',cmap(6,:),'LineWidth',1)
     hold on
     plot(median(distSort(:,1)./distSort(:,2),'omitnan')*[1 1],[0 35],'Color','black','LineWidth',3,'LineStyle','--')
     hold off
    set(gca,'FontSize',24,'LineWidth',1,'FontName','Arial')
    xlabel('Ratio','FontSize',24,'FontName','Arial')
    ylabel('Count','FontSize',24,'FontName','Arial')
    title('Ratio Short/Middle Edges')
    
    % auto-set axes based on maximum bin height + 25%
    bin_counts4 = get(h, 'Values')';
    max_count = max(bin_counts4);
    ylim([0,max_count*1.25])

    % add a text box showing median value
    num = median(distSort(:,1)./distSort(:,2));
    str = sprintf('%.2g', num);
    txt = ['median = ' str];
    annotation('textbox', [0.175, 0.75, 0.1, 0.1], 'String', txt,'FontSize',24,'EdgeColor','none','FontName','Arial');

    box off
    saveas(gcf,[saveNames{kk} '_SmallestDistanceRatioMiddleDistance.pdf'],'pdf')
 %%   

    figure(5)
    h = histogram(distSort(:,1)./distSort(:,3),0:0.1:1.2,'FaceColor',cmap(7,:),'LineWidth',1)
     hold on
     plot(median(distSort(:,1)./distSort(:,3),'omitnan')*[1 1],[0 35],'Color','black','LineWidth',3,'LineStyle','--')
     hold off
    set(gca,'FontSize',24,'LineWidth',1,'FontName','Arial')
    xlabel('Ratio','FontSize',24,'FontName','Arial')
    ylabel('Count','FontSize',24,'FontName','Arial')
    title('Ratio Short/Long Edges')

    % auto-set axes based on maximum bin height + 25%
    bin_counts5 = get(h, 'Values')';
    max_count = max(bin_counts5);
    ylim([0,max_count*1.25])

% add a text box showing median value
    num = median(distSort(:,1)./distSort(:,3));
    str = sprintf('%.2g', num);
    txt = ['median = ' str];
    annotation('textbox', [0.175, 0.75, 0.1, 0.1], 'String', txt,'FontSize',24,'EdgeColor','none','FontName','Arial');

    box off
    saveas(gcf,[saveNames{kk} '_SmallestDistanceRatioLargestDistance.pdf'],'pdf')

   %% 


    figure(51)
    h = histogram(distSort(:,2)./distSort(:,3),0:0.1:1.2,'FaceColor',cmap(8,:),'LineWidth',1)
     hold on
     plot(median(distSort(:,2)./distSort(:,3),'omitnan')*[1 1],[0 35],'Color','black','LineWidth',3,'LineStyle','--')
     hold off
    set(gca,'FontSize',24,'LineWidth',1,'FontName','Arial')
    xlabel('Ratio','FontSize',24,'FontName','Arial')
    ylabel('Count','FontSize',24,'FontName','Arial')
    title('Ratio Middle/Long Edges')

    % auto-set axes based on maximum bin height + 25%
    bin_counts6 = get(h, 'Values')';
    max_count = max(bin_counts6);
    ylim([0,max_count*1.25])

    box off

    % add a text box showing median value
    num = median(distSort(:,2)./distSort(:,3));
    str = sprintf('%.2g', num);
    txt = ['median = ' str];
    annotation('textbox', [0.175, 0.75, 0.1, 0.1], 'String', txt,'FontSize',24,'EdgeColor','none','FontName','Arial');


    saveas(gcf,[saveNames{kk} '_MiddleDistanceRatioLargestDistance.pdf'],'pdf')

    %% Mean Distance

   fig = figure(6)
 h = histogram(mean(distSort,2,'omitnan'),10:5:60,'FaceColor',cmap(1,:),'LineWidth',1)
    set(gca,'FontSize',20,'LineWidth',1,'FontName','Arial')
    xlabel('Distance (nm)','FontSize',20,'FontName','Arial')
    ylabel('Count','FontSize',20,'FontName','Arial')
    title(["Mean", "Inter-Blade Distance"])

    width = fig.Position(3);
    pos = get(gcf, 'Position');
    pos(3) = width/1.33;
    set(gcf, 'Position', pos);

        % auto-set axes based on maximum bin height + 25%
    bin_counts = get(h, 'Values');
    max_count = max(bin_counts);
    ylim([0,max_count*1.25])

    box off
    saveas(gcf,[saveNames{kk} '_MeanDistance.pdf'],'pdf')

%     %% Scatter Plots
% 
%     figure(7)
%     plot(distSort(:,3)-distSort(:,1),distSort(:,3)-distSort(:,2),'o')


end

close all
