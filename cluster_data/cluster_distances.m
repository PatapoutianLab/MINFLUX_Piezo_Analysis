function clusterInfo = cluster_distances(clusterInfo,varargin)

if nargin>1
    plotDistances = varargin{1};
else
    plotDistances = 1;
end

% check which color entries are empty
c_idx = [];
for i=1:3
    if numel(clusterInfo.meanClus)>i-1
        if ~isempty(clusterInfo.meanClus{i})
            c_idx = [c_idx,i];
        end
    end
end


if length(c_idx)==1
    
    [clusterDist.idxDnn1,clusterDist.Dnn1] = knnsearch(clusterInfo.meanClus{c_idx(1)},clusterInfo.meanClus{c_idx(1)},'K',2);
    %     connVec = clusterInfo.meanClus{c_idx(1)}-clusterInfo.meanClus{c_idx(1)}(clusterDist.idxDnn1(:,2),:);
    %     clusterDist.phiDnn1 = atan2(connVec(:,2),connVec(:,1));
    %     clusterDist.thetaDnn1 = acos(connVec(:,3)./arrayfun(norm,connVec));
    clusterDist.distAllC1 = pdist(clusterInfo.meanClus{c_idx(1)});
    %%% Plot distances
    if plotDistances
        edgesNearest = 0:3:300;
        
        fig = figure(445);
        clf
        
        subplot(431)
        %         histogram(clusterDist.phiDnn1)
        xlabel('phi (rad)')
        
        subplot(432)
        %         histogram(clusterDist.thetaDnn1)
        xlabel('theta (rad)')
        
        subplot(437)
        histogram(clusterDist.Dnn1(:,2),edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([median(clusterDist.Dnn1(:,2)),median(clusterDist.Dnn1(:,2))],get(gca,'YLim'),'DisplayName',num2str(median(clusterDist.Dnn1(:,2))))
        xlabel(['color 1 #', num2str(clusterInfo.clusNum{1}),' (nm)'])
        legend show
        
        subplot(438)
        histogram(clusterDist.distAllC1,'EdgeColor','none','FaceColor','b')
        hold on
        plot([median(clusterDist.distAllC1),median(clusterDist.distAllC1)],get(gca,'YLim'),'DisplayName',num2str(median(clusterDist.distAllC1)))
        xlabel('color 1 (nm)')
        legend show
    end
    
    
elseif length(c_idx)==2
    
    %% Calculate distance between two colors
    [~,clusterDist.DnnC1] = knnsearch(clusterInfo.meanClus{c_idx(2)},clusterInfo.meanClus{c_idx(1)});
    [~,clusterDist.DnnC2] = knnsearch(clusterInfo.meanClus{c_idx(1)},clusterInfo.meanClus{c_idx(2)});
    [~,clusterDist.Dnn1] = knnsearch(clusterInfo.meanClus{c_idx(1)},clusterInfo.meanClus{1},'K',2);
    [~,clusterDist.Dnn2] = knnsearch(clusterInfo.meanClus{c_idx(2)},clusterInfo.meanClus{2},'K',2);
    
    clusterDist.crossDistAll = pdist2(clusterInfo.meanClus{c_idx(1)},clusterInfo.meanClus{c_idx(2)});
    
    allClus = [squeeze(clusterInfo.meanClus{c_idx(1)}); clusterInfo.meanClus{c_idx(2)}];
    clusterDist.allClusDist = pdist(allClus);
    
    clusterDist.distAllC1 = pdist(clusterInfo.meanClus{c_idx(1)});
    clusterDist.distAllC2 = pdist(clusterInfo.meanClus{c_idx(2)});
    
    
    %%% Plot distances
    if plotDistances
        edgesNearest = 0:5:300;
        fig = figure(445);
        clf
        
        % nearest neighbour
        subplot(431)
        histogram(clusterDist.DnnC1,edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.DnnC1),nanmedian(clusterDist.DnnC1)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.DnnC1)))
        xlabel(['inter-color 1 #', num2str(clusterInfo.clusNum{c_idx(1)}),' (nm)'])
        legend show
        title('nearest neighbour')
        
        subplot(434)
        histogram(clusterDist.DnnC2,edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.DnnC2),nanmedian(clusterDist.DnnC2)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.DnnC2)))
        xlabel(['inter-color 2 #', num2str(clusterInfo.clusNum{c_idx(2)}),' (nm)'])
        legend show
        % title(['number of clusters ',num2str(length(DnnC2))])
        
        subplot(437)
        histogram(clusterDist.Dnn1(:,2),edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.Dnn1(:,2)),nanmedian(clusterDist.Dnn1(:,2))],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.Dnn1(:,2))))
        xlabel(['color 1 #', num2str(clusterInfo.clusNum{c_idx(1)}),' (nm)'])
        legend show
        
        subplot(4,3,10)
        histogram(clusterDist.Dnn2(:,2),edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.Dnn2(:,2)),nanmedian(clusterDist.Dnn2(:,2))],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.Dnn2(:,2))))
        xlabel(['color 2 #', num2str(clusterInfo.clusNum{c_idx(2)}),' (nm)'])
        legend show
        
        % all
        subplot(432)
        histogram(clusterDist.allClusDist(:),'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.allClusDist(:)),nanmedian(clusterDist.allClusDist(:))],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.allClusDist(:))))
        xlabel('both colors (nm)')
        legend show
        title('all-to-all')
        
        subplot(435)
        histogram(clusterDist.crossDistAll(:),'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.crossDistAll(:)),nanmedian(clusterDist.crossDistAll(:))],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.crossDistAll(:))))
        xlabel('cross colors (nm)')
        legend show
        
        subplot(438)
        histogram(clusterDist.distAllC1,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.distAllC1),nanmedian(clusterDist.distAllC1)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.distAllC1)))
        xlabel('color 1 (nm)')
        legend show
        
        subplot(4,3,11)
        histogram(clusterDist.distAllC2,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.distAllC2),nanmedian(clusterDist.distAllC2)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.distAllC2)))
        xlabel('color 2 (nm)')
        legend show
    end
elseif length(c_idx)==3
    %% Calculate distance between two colors
    [~,clusterDist.DnnC12] = knnsearch(clusterInfo.meanClus{1},clusterInfo.meanClus{2});
    [~,clusterDist.DnnC21] = knnsearch(clusterInfo.meanClus{2},clusterInfo.meanClus{1});
    [~,clusterDist.DnnC13] = knnsearch(clusterInfo.meanClus{1},clusterInfo.meanClus{3});
    [~,clusterDist.DnnC31] = knnsearch(clusterInfo.meanClus{3},clusterInfo.meanClus{1});
    [~,clusterDist.DnnC23] = knnsearch(clusterInfo.meanClus{2},clusterInfo.meanClus{3});
    [~,clusterDist.DnnC32] = knnsearch(clusterInfo.meanClus{3},clusterInfo.meanClus{2});
    [~,clusterDist.Dnn1] = knnsearch(clusterInfo.meanClus{1},clusterInfo.meanClus{1},'K',2);
    [~,clusterDist.Dnn2] = knnsearch(clusterInfo.meanClus{2},clusterInfo.meanClus{2},'K',2);
    [~,clusterDist.Dnn3] = knnsearch(clusterInfo.meanClus{3},clusterInfo.meanClus{3},'K',2);
    
    %     clusterDist.crossDistAll = pdist2(clusterInfo.meanClus{1},clusterInfo.meanClus{2});
    
    allClus = [squeeze(clusterInfo.meanClus{1}); clusterInfo.meanClus{2};squeeze(clusterInfo.meanClus{3})];
    clusterDist.allClusDist = pdist(allClus);
    
    clusterDist.distAllC1 = pdist(clusterInfo.meanClus{1});
    clusterDist.distAllC2 = pdist(clusterInfo.meanClus{2});
    clusterDist.distAllC3 = pdist(clusterInfo.meanClus{3});
    
    %%% Plot distances
    if plotDistances
        edgesNearest = 0:10:300;
        fig = figure(445);
        clf
        
        % nearest neighbour
        subplot(631)
        histogram(clusterDist.DnnC12,edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.DnnC12),nanmedian(clusterDist.DnnC12)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.DnnC12)))
        xlabel(['inter-color 12 #', num2str(clusterInfo.clusNum{2}),' (nm)'])
        legend show
        title('nearest neighbour')
        
        subplot(634)
        histogram(clusterDist.DnnC21,edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.DnnC21),nanmedian(clusterDist.DnnC21)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.DnnC21)))
        xlabel(['inter-color 21 #', num2str(clusterInfo.clusNum{1}),' (nm)'])
        legend show
        % title(['number of clusters ',num2str(length(DnnC2))])
        
        subplot(637)
        histogram(clusterDist.DnnC13,edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.DnnC13),nanmedian(clusterDist.DnnC13)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.DnnC13)))
        xlabel(['inter-color 13 #', num2str(clusterInfo.clusNum{3}),' (nm)'])
        legend show
        title('nearest neighbour')
        
        subplot(6,3,10)
        histogram(clusterDist.DnnC31,edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.DnnC31),nanmedian(clusterDist.DnnC31)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.DnnC31)))
        xlabel(['inter-color 31 #', num2str(clusterInfo.clusNum{1}),' (nm)'])
        legend show
        % title(['number of clusters ',num2str(length(DnnC2))])
        
        subplot(6,3,13)
        histogram(clusterDist.DnnC23,edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.DnnC23),nanmedian(clusterDist.DnnC23)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.DnnC23)))
        xlabel(['inter-color 23 #', num2str(clusterInfo.clusNum{3}),' (nm)'])
        legend show
        title('nearest neighbour')
        
        subplot(6,3,16)
        histogram(clusterDist.DnnC32,edgesNearest,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.DnnC32),nanmedian(clusterDist.DnnC32)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.DnnC32)))
        xlabel(['inter-color 32 #', num2str(clusterInfo.clusNum{2}),' (nm)'])
        legend show
        % title(['number of clusters ',num2str(length(DnnC2))])
        
        try
            subplot(6,3,2)
            histogram(clusterDist.Dnn1(:,2),edgesNearest,'EdgeColor','none','FaceColor','b')
            hold on
            plot([nanmedian(clusterDist.Dnn1(:,2)),nanmedian(clusterDist.Dnn1(:,2))],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.Dnn1(:,2))))
            xlabel(['color 1 #', num2str(clusterInfo.clusNum{1}),' (nm)'])
            legend show
        catch
        end
        
        try
            subplot(6,3,5)
            histogram(clusterDist.Dnn2(:,2),edgesNearest,'EdgeColor','none','FaceColor','b')
            hold on
            plot([nanmedian(clusterDist.Dnn2(:,2)),nanmedian(clusterDist.Dnn2(:,2))],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.Dnn2(:,2))))
            xlabel(['color 2 #', num2str(clusterInfo.clusNum{2}),' (nm)'])
            legend show
        catch
        end
        
        try
            subplot(6,3,8)
            histogram(clusterDist.Dnn3(:,2),edgesNearest,'EdgeColor','none','FaceColor','b')
            hold on
            plot([nanmedian(clusterDist.Dnn3(:,2)),nanmedian(clusterDist.Dnn3(:,2))],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.Dnn3(:,2))))
            xlabel(['color 2 #', num2str(clusterInfo.clusNum{3}),' (nm)'])
            legend show
        catch
        end
        
        
        %%
        % all
        subplot(633)
        histogram(clusterDist.allClusDist(:),'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.allClusDist(:)),nanmedian(clusterDist.allClusDist(:))],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.allClusDist(:))))
        xlabel('both colors (nm)')
        legend show
        title('all-to-all')
        
        %     subplot(636)
        %     histogram(clusterDist.crossDistAll(:),'EdgeColor','none','FaceColor','b')
        %     hold on
        %     plot([nanmedian(clusterDist.crossDistAll(:)),nanmedian(clusterDist.crossDistAll(:))],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.crossDistAll(:))))
        %     xlabel('cross colors (nm)')
        %     legend show
        
        subplot(639)
        histogram(clusterDist.distAllC1,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.distAllC1),nanmedian(clusterDist.distAllC1)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.distAllC1)))
        xlabel('color 1 (nm)')
        legend show
        
        subplot(6,3,12)
        histogram(clusterDist.distAllC2,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.distAllC2),nanmedian(clusterDist.distAllC2)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.distAllC2)))
        xlabel('color 2 (nm)')
        legend show
        
        subplot(6,3,15)
        histogram(clusterDist.distAllC3,'EdgeColor','none','FaceColor','b')
        hold on
        plot([nanmedian(clusterDist.distAllC3),nanmedian(clusterDist.distAllC3)],get(gca,'YLim'),'DisplayName',num2str(nanmedian(clusterDist.distAllC3)))
        xlabel('color 3 (nm)')
        legend show
    end
else
    disp('Too many colors!')
end


%% Plot combined scatter plot

colordata{1} = 'r';
colordata{2} = 'g';
colordata{3} = 'b';

colorcross{1} = 'r';
colorcross{2} = 'g';
colorcross{3} = 'b';

legendName{1} = 'color1';
legendName{2} = 'color2';
legendName{3} = 'color3';

if plotDistances
    fig = figure(446);
    clf
    for ll=1:numel(clusterInfo.meanClus)
        if clusterInfo.prmAll.cluster.dim==3
            scatter3(clusterInfo.locAll{ll}(:,1),clusterInfo.locAll{ll}(:,2),clusterInfo.locAll{ll}(:,3),30,colordata{ll},'.','DisplayName',legendName{ll})
            hold on
            scatter3(clusterInfo.meanClus{ll}(:,1),clusterInfo.meanClus{ll}(:,2),clusterInfo.meanClus{ll}(:,3),70,colorcross{ll},'x','LineWidth',3,'DisplayName',legendName{ll})
        else
            scatter(clusterInfo.locAll{ll}(:,1),clusterInfo.locAll{ll}(:,2),30,colordata{ll},'.','DisplayName',legendName{ll})
            scatter(clusterInfo.meanClus{ll}(:,1),clusterInfo.meanClus{ll}(:,2),70,colorcross{ll},'x','LineWidth',3,'DisplayName',legendName{ll})
        end
        axis equal
        legend show
        xlabel('x (nm)')
        ylabel('y (nm)')
        if size(clusterInfo.locAll{ll},2)==3
            zlabel('z (nm)')
        end
    end
    clusterInfo.clusterDist = clusterDist;
end


