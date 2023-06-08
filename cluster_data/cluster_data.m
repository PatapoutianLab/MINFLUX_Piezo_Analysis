function clusterInfo = cluster_data(locAll,varargin)

% the returned index contains noise data with idx=1
% locAll %um
prmAll = varargin{2};
prmAll.cluster = set_default_parameters(varargin{:});
prm = prmAll.cluster;

for ll=1:numel(locAll)
    if ~isempty(locAll{ll})
        locAll{ll} = locAll{ll}*1e3;
        loc = locAll{ll}(:,1:prm.dim);
        %% Find clusters
        switch prm.algorithm
            case 'meanShift'
                [~,idx,~] = MeanShiftCluster(loc,prm.meanShift.bandWidth,1);
            case 'dbscan'  % density-based clustering
                [idx,isnoise] = DBSCAN(loc,prm.dbscan.eps,prm.dbscan.minPts);
            case 'dbscan2'
                [idx] = dbscan2(loc, 'params', prm.dbscan2);
            case 'dbscanCov'
                [idx] = dbscanCov(loc, 'params', prm.dbscanCov);
            case 'kmeans' % Non-hierarchical clustering (k-means)
                idx = kmeans(loc,prm.kmeans.numClus);
            case 'gaussianMixture' % Gaussian mixture model %
                [idx] = kmeans(loc,prm.gmm.numClus);
                options = statset('MaxIter',1000);
                gmfit = fitgmdist(loc,prm.gmm.numClus,'Options',options,'Start',idx);
                idx = cluster(gmfit,loc);
            case 'hierarch'
                %                 idx = clusterdata(loc,prm.hierarch.cutoff);
                Y = pdist(loc, 'euclid');
                Z = linkage(Y, 'average');
                %                 inc = inconsistent(Z);
                %                 c = cophenet(Z,Y);
                %                 disp(num2str((inc(:,4))))
                idx = cluster(Z, 'cutoff', prm.hierarch.cutoff);
            case 'spectral'
                idx = spectralcluster(loc,prm.spectral.numClus);
                
            case 'iter'
                meanSigma = prm.dbscan.eps;
                meanNump = prm.dbscan.minPts;
                clear msig
                clear numclus
                for kk=1:50
                    idx = DBSCAN(loc,meanSigma,meanNump);
                    [idx,C,sumd] = kmeans(loc,max(idx+1));
                    for jj=1:max(idx)
                        sigma(jj) = mean(std(loc(idx==jj)));
                        nump(jj) = sum(idx==jj);
                    end
                    meanSigma = mean(sigma);
                    msig(kk) = meanSigma;
                    numclus(kk) = max(idx);
                    meanNump = round(mean(nump));
                    disp(['sigma ',num2str(meanSigma)])
                    disp(['nump ', num2str(meanNump)])
                    disp(['numClus ',num2str(max(idx))])
                    figure(123)
                    clf
                    subplot(211)
                    plot(msig)
                    ylabel('sigma_{clus} (nm)')
                    subplot(212)
                    plot(numclus)
                    xlabel('iterations')
                    ylabel('# loc/cluster')
                end
                
        end
%         disp(['%% ',prm.algorithm, ' %%'])
%         disp(['Number of clusters: ',num2str(max(idx))])
        
        %% Extract cluster information
        
        numClus = max(idx)-min(idx)+1;
        if min(idx)==0
            idx = idx+1;
        elseif min(idx)==-1
            idx = idx+2;
        end
        for jj=min(idx):max(idx)
            if prm.filterNoise
                locClus{ll,jj} = loc((idx==jj)&~isnoise,:);
                meanClus{ll}(jj,:) = mean(loc((idx==jj)&~isnoise,:)',2);
                stdClus{ll}(jj,:) = std(loc((idx==jj)&~isnoise,:)',0,2);
                locAllNoise{ll} = loc(isnoise,:);
                numEl{ll}(jj) = sum(idx==jj);
                
            else
                locClus{ll,jj} = loc(idx==jj,:);
                meanClus{ll}(jj,:) = mean(loc(idx==jj,:)',2);
                stdClus{ll}(jj,:) = std(loc(idx==jj,:)',0,2);
                locAllNoise{ll} = loc;
                numEl{ll}(jj) = sum(idx==jj);
            end
        end
        distAll{ll} = pdist(meanClus{ll});
        clusIdx{ll} = idx;
        clusNum{ll} = numClus;
        clusIso{ll} = min(stdClus{ll},[],2)./max(stdClus{ll},[],2);

    end
end
clusterInfo.locAll = locAll;
clusterInfo.locClus = locClus;
clusterInfo.meanClus = meanClus;
clusterInfo.prmAll = prmAll;
clusterInfo.stdClus = stdClus;
clusterInfo.locAllNoise = locAllNoise;
clusterInfo.clusIdx = clusIdx;
clusterInfo.clusNum = clusNum;
clusterInfo.numEl = numEl;
clusterInfo.clusIso = clusIso;
if prm.plot
    plot_clusters(clusterInfo);
end
end

function pout = set_default_parameters(varargin)

pout = struct;

%         % Data selection
%         %%%%%%%%%%%%%%%%%%%%%%%%%%
%         p.colors = 'all'; % 'color1','color2','color3','all','combined'
%         p.data = 'all'; % 'selection', 'all'

pout.dim = 3;
pout.algorithm = 'dbscan';
pout.filterNoise = 0;
pout.plot = 1;

% Algorithm parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dbscan
pout.dbscan.eps = 15;
pout.dbscan.minPts = 5;
pout.dbscan.distance = 'euclidean';
% dbscan2
pout.dbscan2.db1_size = 20;
pout.dbscan2.db1_minPts = 4;
pout.dbscan2.db2_size = 8;
pout.dbscan2.db2_minPts = 4;
pout.dbscan2.gmm_sigma = 5;
pout.dbscan2.quiet = 1;
% meanShift
pout.meanShift.bandWidth = 10;
% kmeans
pout.kmeans.numClus = 80;
% hierarch
pout.hierarch.cutoff = 1.1547;
% Gaussian mixture
pout.gmm.numClus = 80;
% spectral clustering
pout.spectral.numClus = 80;

%##########################################################################
% Override default parameters with name/value pairs from the input.
if numel(varargin) > 0
    pin = varargin{2}.cluster;
    pout = omex_read_params(pout, {'params',pin});
end

end
