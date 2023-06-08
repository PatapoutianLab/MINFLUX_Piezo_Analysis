
% Code to analyze triple-labeled PIEZO1 molecules from Abberior Instruments 3D MINFLUX localizations
% Written by Eric Mulhall
% cluster_data.m and dependencies written by Pape, J.K., et al: https://doi.org/10.1073/pnas.2009364117

% Instructions:
% First, export valid final localizations from the MINFLUX Imspector interface as a .mat file
% Import the .mat file into the current path
% If necessary, determine proper Z-filtering thresholds from examining localizations in Paraview or Matlab

% Results:
% The final identified molecules will be in the output table 'Results'

%Path to files
%Path = '/Users/path/';

% before starting, add cluster_data.m and dependencies to path

close all
clear
clc

%% Load localization data 

% load .mat file
mfx_data = ' .mat'
load(mfx_data);
% set number format (optional)
%format shortE

%% Filter localization data

% isolate localization and parameter arrays

% localizations
itr_loc = itr.loc;
% counts at offset
itr_eco = itr.eco;
% counts at center position
itr_ecc = itr.ecc;
% center frequency ratio
itr_cfr = itr.cfr;
% frequency at offset
itr_efo = itr.efo;

% define x,y,z coordinates for each localization and link the trace ID
x = itr_loc(:,10,1);
y = itr_loc(:,10,2);
z = itr_loc(:,10,3);
xyz = [x y z];
% Convert tid into double for import
tid = double(tid)';
traces = [x(:,1),y(:,1),z(:,1),tid(:,1)];

% input Z-filt threshold if required
%z_threshold = 0e-07

% Set standard deviation threshold per trace for each dimension (default = 10 nm)
stdev_trace_threshold = 10e-09

% Set # localizations per trace threshold (defalut > 3)
loc_per_trace_threshold = 3

% filter #1: Z threshold
traces_filt = traces;
% output traces greater than Z-threshold
%traces_filt(traces_filt(:, 3) < z_threshold, :)= [];
% output traces less than Z-threshold
%traces_filt(traces_filt(:, 3) > z_threshold, :)= [];

% Filter #2: standard deviation per trace in each dimension
[uv_tid, ~, id_tid] = unique(traces_filt(:,4));
stdev_trace_x = [accumarray(id_tid,traces_filt(:,1),[],@std)];
traces_filt = traces_filt(ismember(traces_filt(:,4), uv_tid(stdev_trace_x < stdev_trace_threshold)),:);
% 
[uv_tid, ~, id_tid] = unique(traces_filt(:,4));
stdev_trace_y = [accumarray(id_tid,traces_filt(:,2),[],@std)];
traces_filt = traces_filt(ismember(traces_filt(:,4), uv_tid(stdev_trace_y < stdev_trace_threshold)),:);
% 
[uv_tid, ~, id_tid] = unique(traces_filt(:,4));
stdev_trace_z = [accumarray(id_tid,traces_filt(:,3),[],@std)];
traces_filt = traces_filt(ismember(traces_filt(:,4), uv_tid(stdev_trace_z < stdev_trace_threshold)),:);
 
% Filter #3: localizations per trace
% Unique elements and locations in third column
[uv_tid, ~, id_tid] = unique(traces(:,4));
% How many of each?
n_tid = histcounts(id_tid,"BinWidth",1);
% Keep ones with more than 5.
traces_filt = traces_filt(ismember(traces_filt(:,4), uv_tid(n_tid > loc_per_trace_threshold)),:);

% Plot figure to check that Z-filtering is correct
figure
scatter3(traces_filt(:,1),traces_filt(:,2),traces_filt(:,3),5,traces_filt(:,4))
title('z-filtering check')

% Convert units for compatibility with later graphs
traces_filt_nano = traces_filt*10^9;

% plot the filtered data in nanometer-scale colored by tid
figure
scatter3(traces_filt_nano(:,1),traces_filt_nano(:,2),traces_filt_nano(:,3),5,traces_filt_nano(:,4))
title('filtered traces nanometer scale')
axis equal
hold on

% Put XYZ coordinates into filt_XYZ for cluster analysis 
filt_XYZ = traces_filt(:,1:3);

% Put filt_XYZ into locAll for cluster analysis
locAll = {[filt_XYZ]*10^6};


%% Define dbscan and GMM parameters and load the cluster algorithm

p.cluster.algorithm = 'dbscan2';
p.cluster.dbscan2.db1_size = 30; 
p.cluster.dbscan2.db1_minPts = 5;
% db2_size is the only parameter changed, between 6-7 depending on noise
p.cluster.dbscan2.db2_size = 7;
p.cluster.dbscan2.db2_minPts = 5;
p.cluster.dbscan2.gmm_sigma = 5;
p.cluster.plot = 1;

p_save.plot_save=1;
p_save.plot_save_path = '04_Plots';
p_save.plot_save_filename_start = '';%p.cluster.algorithm;


%% Perform the 2-step clustering algorithm and GMM fitting

% cluster the data
clusInfo{1} = cluster_data(locAll,'params',p);
        % note that the cluster_data algorithm does not like the array size too large. In that case, the data will need to be split up between locAll{1...X}
        prmPlot.coloring = 'NN';
        clusInfo{1} = cluster_distances(clusInfo{1},0);

clust_XYZ = clusInfo{1, 1}.meanClus{1, 1};
size_clust = size(clust_XYZ);
size_clust = size_clust(:,1);
clust_index = [];
clust_index = (1:size_clust)';
clust_XYZ = [clust_XYZ clust_index];
clust_stdev = clusInfo{1, 1}.stdClus{1, 1};
% average stdev in 3 dimensions for plotting
clust_stdev_avg = ([(clust_stdev(:,1) + clust_stdev(:,2) + clust_stdev(:,3))/3]);

% create an XYZ + average stdev array for visualization
clust_XYZ_stdev_avg = [clust_XYZ clust_stdev_avg];

% create an XYZ + XYZ_stdev array
clust_XYZ_stdev = [clust_XYZ clust_stdev];

% plot the filtered data and the clustered data
figure
% plot cluster center positions
scatter3(clust_XYZ_stdev_avg(:,1),clust_XYZ_stdev_avg(:,2),clust_XYZ_stdev_avg(:,3),100,"black")
axis equal
hold on
% plot filtered traces
scatter3(traces_filt_nano(:,1),traces_filt_nano(:,2),traces_filt_nano(:,3),5,traces_filt_nano(:,4))
colormap jet
title('Identified Clusters')

%% DBSCAN the clusters to look for clusters of 3

% define the epsilon for dbscan (default = 100)
epsilon = 100;
% define minpts for dbscan (default = 3)
minpts = 3;

% % optional: estimate epsilon
% kD = pdist2(clust_XYZ,clust_XYZ,'euc','Smallest',minpts);
% figure
% plot(sort(kD(end,:)));
% title('k-distance graph')
% xlabel('Points sorted with 50th nearest distances')
% ylabel('50th nearest distances')
% grid

% run the DBSCAN for XYZ coordinates
idx = dbscan(clust_XYZ_stdev(:,1:3),epsilon,minpts);

% create a merged array with XYZ, StdDev, and DBSCAN results
clust_XYZ_DBSCAN = [clust_XYZ_stdev idx];

% remove rows that have a -1 value, i.e. they are not in a cluster with minimum 3 points within epsilon
clust_XYZ_DBSCAN(clust_XYZ_DBSCAN(:, 8)== -1, :)= [];

% plot the data to check for proper removal of points
figure
scatter3(clust_XYZ_DBSCAN(:,1),clust_XYZ_DBSCAN(:,2),clust_XYZ_DBSCAN(:,3),"black","filled","MarkerFaceAlpha",0.6)
axis equal
hold on
title ('raw DBSCAN results')

%% Eliminate non-trimers from DBSCAN clusters

% Unique elements and locations in third column
[uv, ~, id] = unique(clust_XYZ_DBSCAN(:,8));
% How many of each?
n = histcounts(id,'BinWidth',1);
% Keep traces where cluster number = 3
clust_XYZ_DBSCAN_3 = clust_XYZ_DBSCAN(ismember(clust_XYZ_DBSCAN(:,8), uv(n==3)),:);

% sort the rows
clust_XYZ_DBSCAN_3_sort = sortrows(clust_XYZ_DBSCAN_3,8);

% plot clusters for sanity check
figure
scatter3(clust_XYZ_DBSCAN_3_sort(:,1),clust_XYZ_DBSCAN_3_sort(:,2),clust_XYZ_DBSCAN_3_sort(:,3),100,'black','LineWidth',1)
axis equal
hold on
% plot filtered traces
scatter3(traces_filt_nano(:,1),traces_filt_nano(:,2),traces_filt_nano(:,3),5,traces_filt_nano(:,4))
colormap jet
title('DBSCAN trimers identified')

%% Nearest neighbor analysis of all data

% this is a separate nearest neighbor analysis of all data that will then be compared to the DBSCAN results

% define nearest neighbor parameters
neighborNumber = 2; % Number of neighbors each peak should have
neighborDist = 50; % nm, maximum neighbor distance
minDist = 6; % nm; minimum neighbor distance. Zero gets rid of self comparisons

Z_dist = squareform(pdist(clust_XYZ_stdev(:,1:3))); % all distances between neighbors
Zplus2 = (Z_dist < neighborDist) & (Z_dist > minDist);

% set a requirement on neighbor number, export a logical, then turn into double
plus2 = sum(Zplus2,2) == neighborNumber;
plus2 = double(plus2);

% create a merged array with XYZ, average StdDev, and NN results
clust_XYZ_NN = [clust_XYZ_stdev plus2];

% remove rows that have a 0 value, i.e. they fail the NN requirements
clust_XYZ_NN(clust_XYZ_NN(:, 8)== 0, :)= [];

figure
scatter3(clust_XYZ_NN(:,1),clust_XYZ_NN(:,2),clust_XYZ_NN(:,3),100,'black','LineWidth',1)
axis equal
hold on
% plot filtered traces
scatter3(traces_filt_nano(:,1),traces_filt_nano(:,2),traces_filt_nano(:,3),5,traces_filt_nano(:,4))
colormap jet
title('Clusters meeting NN requirements')

%% Match the NN and DBSCAN data based on trace ids

% compare cluster IDs between the DBSCAN and Nearest Neighbor analyses
[~, ind] = ismember(clust_XYZ_DBSCAN_3_sort(:,4), clust_XYZ_NN(:,4),'rows');
clust_comp_DBSCAN_NN = [clust_XYZ_DBSCAN_3_sort ind];
% remove clusters not meeting requirements from both analyses
clust_comp_DBSCAN_NN(clust_comp_DBSCAN_NN(:, 9)== 0, :)= [];

% make sure that the candidate trimers contain only 3 clusters
% Count the number of occurrences of each DBSCAN ID
counts = accumarray(clust_comp_DBSCAN_NN(:,8),1);

% Find the numbers that have 3 occurrences
numbersWith3occurrences = find(counts==3);
% remove those with fewer than 3 occurrences
clust_comp_DBSCAN_NN_3 = clust_comp_DBSCAN_NN(ismember(clust_comp_DBSCAN_NN(:,8),numbersWith3occurrences(:,1)),:);

% Rearrange array in this order: X, Y, Z, Xstdev, Ystdev, Zstdev, DBSCAN_ID
clust_comp_DBSCAN_NN_3_sort = clust_comp_DBSCAN_NN_3(:,[1:3 5:7 8]);

% sort rows based on cluster ID
clust_comp_DBSCAN_NN_3_sort = sortrows(clust_comp_DBSCAN_NN_3_sort,7);

% plot final identified trimers meeting all requirements over the filtered traces
figure
scatter3(traces_filt_nano(:,1),traces_filt_nano(:,2),traces_filt_nano(:,3),5,traces_filt_nano(:,4))
hold on
scatter3(clust_comp_DBSCAN_NN_3_sort(:,1),clust_comp_DBSCAN_NN_3_sort(:,2),clust_comp_DBSCAN_NN_3_sort(:,3),50,"red","filled","MarkerFaceAlpha",0.6)
axis equal
title('Identified trimers meeting NN and DBSCAN requirements')

%% Calculate angles between each trimer 

% loop over each grouping of three positions (trimer)
numtrimers = size(clust_comp_DBSCAN_NN_3_sort, 1) / 3;
for i = 1:numtrimers
    % calculate the starting index for the current trimer
    startIndex = (i - 1) * 3 + 1;
    % extract the positions for the current trimer
    cluster = clust_comp_DBSCAN_NN_3_sort(startIndex:startIndex+2, :);
    
    % loop over each point in the timer
    for j = 1:3
        % starting index for the current point
        pointIndex = startIndex + j - 1;
        
        % indices for the other two points in the trimer
        otherPointIndex1 = startIndex + mod(j, 3);
        otherPointIndex2 = startIndex + mod(j + 1, 3);
        
        vector1 = cluster(otherPointIndex1 - startIndex + 1, :) - cluster(pointIndex - startIndex + 1, :);
        vector2 = cluster(otherPointIndex2 - startIndex + 1, :) - cluster(pointIndex - startIndex + 1, :);
        
        % normalize and calculate angles
        normalizedVector1 = vector1 / norm(vector1);
        normalizedVector2 = vector2 / norm(vector2);
        dotProduct = dot(normalizedVector1, normalizedVector2);
        angleRadians = acos(dotProduct);
        angleDegrees = rad2deg(angleRadians);
        
        % add the angle to the array of trimers
        clust_comp_DBSCAN_NN_3_sort(pointIndex, 8) = angleDegrees;
    end
end

%% Remove trimers that have angles > 120 degrees

clusterIDs = unique(clust_comp_DBSCAN_NN_3_sort(:, 7));
clustersToRemove = false(size(clusterIDs));

for i = 1:numel(clusterIDs)
    clusterIndices = (clust_comp_DBSCAN_NN_3_sort(:, 7) == clusterIDs(i));
        if any(clust_comp_DBSCAN_NN_3_sort(clusterIndices, 8) > 120)
        clustersToRemove(i) = true;
    end
end

% remove the trimers with angles >120 degrees
clust_comp_DBSCAN_NN_3_sort_120_filtered = clust_comp_DBSCAN_NN_3_sort(~ismember(clust_comp_DBSCAN_NN_3_sort(:, 7), clusterIDs(clustersToRemove)), :);


%% Calculate inter-blade distances

clusterIDs = unique(clust_comp_DBSCAN_NN_3_sort_120_filtered(:, 7));

for i = 1:numel(clusterIDs)
    % indices of points for trimer
    clusterIndices = (clust_comp_DBSCAN_NN_3_sort_120_filtered(:, 7) == clusterIDs(i));
    
    % XYZ positions of the points
    clusterPoints = clust_comp_DBSCAN_NN_3_sort_120_filtered(clusterIndices, 1:3);
    
    % calculate pairwise distance
    numPoints = size(clusterPoints, 1);
    
    % distances for each point in the cluster
    for k = 1:numPoints
        pointIndices = mod((k-1:k) - 1, numPoints) + 1;
        pointPair = clusterPoints(pointIndices, :);
        pointDistances = sqrt(sum(diff(pointPair).^2, 2));
        assignIndices = find(clusterIndices);
        startIdx = assignIndices(k);
        
        % assign the distances
        clust_comp_DBSCAN_NN_3_sort_120_filtered(startIdx:startIdx, 9) = pointDistances;
    end
end

%% Plot the final identified trimeric PIEZO1 molecules meeting all requirments

figure
scatter3(traces_filt_nano(:,1),traces_filt_nano(:,2),traces_filt_nano(:,3),5,traces_filt_nano(:,4))
hold on
scatter3(clust_comp_DBSCAN_NN_3_sort_120_filtered(:,1),clust_comp_DBSCAN_NN_3_sort_120_filtered(:,2),clust_comp_DBSCAN_NN_3_sort_120_filtered(:,3),50,"red","filled","MarkerFaceAlpha",0.6)
axis equal
title('Identified trimers meeting all requirements')


%% Create table of final results

headers = {'X', 'Y', 'Z', 'stdev X', 'stdev Y', 'stdev Z', 'DBSCAN ID', 'Inter-blade angle', 'Inter-blade distance'};
Results = array2table(clust_comp_DBSCAN_NN_3_sort_120_filtered, 'VariableNames', headers);

%% Save work

save_name = mfx_data;
save_name = save_name(1:end-4);

save([save_name '_analyzed.mat'])

disp('Processing completed!')
