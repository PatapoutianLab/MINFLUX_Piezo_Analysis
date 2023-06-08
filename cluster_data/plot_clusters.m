function [fig_clust] = plot_clusters(clusterInfo,varargin)


p = set_default_parameters(varargin{:});
for ll=1:numel(clusterInfo.locAll)
    if ~isempty(clusterInfo.locAll{ll})
        edges = [0:0.5:20 Inf];
        edgesNumEl = [0:1:100 Inf];
        
        
        fig_clust(ll) = figure(345+ll);
        clf
        
        subplot(3,5,1:10) % localization data
        switch p.coloring
            case 'distinguishable'
                colors = distinguishable_colors(clusterInfo.clusNum{ll});
                tmp = clusterInfo.clusIdx{ll};
                tmp(isnan(tmp)) = 1;
                clusColors = colors(tmp,:);
                if clusterInfo.prmAll.cluster.dim==3
                    scatter3(clusterInfo.locAll{ll}(:,1),clusterInfo.locAll{ll}(:,2),clusterInfo.locAll{ll}(:,3),50,clusColors,'.'), hold on
                    scatter3(clusterInfo.meanClus{ll}(:,1),clusterInfo.meanClus{ll}(:,2),clusterInfo.meanClus{ll}(:,3),50,colors,'x','LineWidth',3)
                else
                    scatter(clusterInfo.locAll{ll}(:,1),clusterInfo.locAll{ll}(:,2),50,clusColors,'.'), hold on
                    scatter(clusterInfo.meanClus{ll}(:,1),clusterInfo.meanClus{ll}(:,2),50,colors,'x','LineWidth',3)
                end
            case 'NN' % use average distance to first nearest neigbours within the same color
                [~,dNN] = knnsearch(clusterInfo.meanClus{ll},clusterInfo.meanClus{ll},'K',3);
                dNN = mean(dNN(:,2:3),2);
                if clusterInfo.prmAll.cluster.dim==3
                    scatter3(clusterInfo.locAll{ll}(:,1),clusterInfo.locAll{ll}(:,2),clusterInfo.locAll{ll}(:,3),50,dNN(clusterInfo.clusIdx{ll}),'.'), hold on
                    scatter3(clusterInfo.meanClus{ll}(:,1),clusterInfo.meanClus{ll}(:,2),clusterInfo.meanClus{ll}(:,3),50,dNN,'x','LineWidth',3);
                else
                    scatter(clusterInfo.locAll{ll}(:,1),clusterInfo.locAll{ll}(:,2),50,dNN(clusterInfo.clusIdx{ll}),'.'), hold on
                    scatter(clusterInfo.meanClus{ll}(:,1),clusterInfo.meanClus{ll}(:,2),50,dNN,'x','LineWidth',3);
                end
                colormap('jet')
            case 'NNcolor' % use nearest neighbour distance to other color
                if ll==1
                    [~,dNN] = knnsearch(clusterInfo.meanClus{2},clusterInfo.meanClus{1});
                elseif ll==2
                    [~,dNN] = knnsearch(clusterInfo.meanClus{1},clusterInfo.meanClus{2});
                end
                if clusterInfo.prmAll.cluster.dim==3
                    scatter3(clusterInfo.locAll{ll}(:,1),clusterInfo.locAll{ll}(:,2),clusterInfo.locAll{ll}(:,3),50,dNN(clusterInfo.clusIdx{ll}),'.'), hold on
                    scatter3(clusterInfo.meanClus{ll}(:,1),clusterInfo.meanClus{ll}(:,2),clusterInfo.meanClus{ll}(:,3),50,dNN,'x','LineWidth',3);
                else
                    scatter(clusterInfo.locAll{ll}(:,1),clusterInfo.locAll{ll}(:,2),50,dNN(clusterInfo.clusIdx{ll}),'.'), hold on
                    scatter(clusterInfo.meanClus{ll}(:,1),clusterInfo.meanClus{ll}(:,2),50,dNN,'x','LineWidth',3);
                end
                colormap('parula')
        end
        
        xlabel('x (nm)')
        ylabel('y (nm)')
        if clusterInfo.prmAll.cluster.dim==3
            zlabel('z (nm)')
        end
        title(['numClus ',num2str(clusterInfo.clusNum{ll})])
        axis equal
        hold on
        
        
        subplot(3,5,11) % isotropy
        %         isotr = min(clusterInfo.stdClus{ll},[],2)./max(clusterInfo.stdClus{ll},[],2);
        edges_isotr = 0:0.05:1;
        histogram(clusterInfo.clusIso{ll},edges_isotr,'EdgeColor','none','FaceColor','b')
        hold on
        med = median(clusterInfo.clusIso{ll});
        plot([med,med],get(gca,'YLim'),'k--','LineWidth',2)
        title(['med ',num2str(med,'%.2f')])
        xlabel('isotr')
        
        subplot(3,5,12) % stdx
        histogram(clusterInfo.stdClus{ll}(:,1),edges,'EdgeColor','none','FaceColor','b')
        hold on
        med = median(clusterInfo.stdClus{ll}(:,1));
        plot([med,med],get(gca,'YLim'),'k--','LineWidth',2)
        title(['med ',num2str(med,'%.1f'), 'nm'])
        xlabel('std_x (nm)')
        
        
        subplot(3,5,13) % stdy
        histogram(clusterInfo.stdClus{ll}(:,2),edges,'EdgeColor','none','FaceColor','b')
        hold on
        med = median(clusterInfo.stdClus{ll}(:,2));
        plot([med,med],get(gca,'YLim'),'k--','LineWidth',2)
        title(['med ',num2str(med,'%.1f'), 'nm'])
        xlabel('std_y (nm)')
        
        if clusterInfo.prmAll.cluster.dim==3
            subplot(3,5,14) % stdz
            histogram(clusterInfo.stdClus{ll}(:,3),edges,'EdgeColor','none','FaceColor','b')
            hold on
            med = median(clusterInfo.stdClus{ll}(:,3));
            plot([med,med],get(gca,'YLim'),'k--','LineWidth',2)
            title(['med ',num2str(med,'%.1f'), 'nm'])
            xlabel('std_z (nm)')
        end
        
        subplot(3,5,15) 
        histogram(clusterInfo.numEl{ll},edgesNumEl,'EdgeColor','none','FaceColor','b')
        hold on
        med = median(clusterInfo.numEl{ll});
        plot([med,med],get(gca,'YLim'),'k--','LineWidth',2)
        title(['med ',num2str(med,'%.0f')])
        xlabel('#loc in cluster (nm)')
    end
    
end

end


function p = set_default_parameters(varargin)

p = struct;
p.coloring = 'distinguishable';
%##########################################################################
% Override default parameters with name/value pairs from the input.
if numel(varargin) > 0
    p = omex_read_params(p, varargin);
end
end

