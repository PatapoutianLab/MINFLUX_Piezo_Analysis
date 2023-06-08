function [c2_idx] = dbscan2(x, varargin)
% dbscan2 - Clustering by two sequential DBSCANs. The first one is applied
% to the points in 'x'. The second one is applied to each cluster
% recognized by the first DBSCAN. 
% 
% After the second DBSCAN, a gaussian mixture classification is applied.
% The mean of each population is obtained from the second DBSCAN clustering
% results. The covariance is a fixed parameter.
% 
% Syntax: 
%  [c2_idx] = function_name(x)
%  function_name(x, 'params', p)
%  function_name(x, Name, Val, ...)
%
% Inputs:
%    x - Points of dimension k to cluster. Array of n x k.
%
% Outputs:
%    c2_idx - Cluster id per point.
%
% Properties
%         p.db1_size = 20;
%         p.db1_minPts = 4;
% 
%         p.db2_size = 7;
%         p.db2_minPts = 4;
% 
%         p.gmm_sigma = 3;
% 
%         p.ground_truth = [];
% 
%         p.quiet = 1;
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: DBSCAN,  OTHER_FUNCTION_NAME2

% Author: 
% email address: 
% Website: 
% December 1999; 
% Last revision: 12-May-2004

p = set_default_parameters(varargin{:});

%##########################################################################
% CLUSTERING 1

% Calculate first DBSCAN clustering
c1_idx = DBSCAN(x, p.db1_size, p.db1_minPts);


% Initialize
c=0;
c2_idx = zeros(size(c1_idx));
c_sub_idx = zeros(size(c1_idx));  % Subclustering results after GMM. Contains 'local indices'

% Run through the clusters
for i=1:max(c1_idx)
    if ~p.quiet, display_progress(i, max(c1_idx), 'waitbar', 10), end

    % Elements in current cluster
    ind_aux = c1_idx==i;
    x_aux = x(ind_aux,:);
    
    % Calculate second DBSCAN clustering
    cl_idx_aux = DBSCAN(x_aux, p.db2_size, p.db2_minPts);
   
    % If is not empty, recluster with a gasusian mixture to not lose
    % points.
    if any(cl_idx_aux)
        gmm_idx = gmm_eval(x_aux, cl_idx_aux, p.gmm_sigma);
        c_sub_idx(ind_aux) = gmm_idx;
        
        
        % Renumber all clusters with a floating index c
        for j=1:max(gmm_idx)
            c = c+1;
            c2_idx(ind_aux & c_sub_idx==j) = c;
        end
    else
        % If the second DBSCAN didn't find it because of low density, I
        % keep it anyway.
        c = c+1;
        c2_idx(ind_aux) = c;
    end
end

end


function [gmm_idx] = gmm_eval(x, idx, sigma)
% Evaluate gaussian mixture model on points X. 
% idx contains the previous clustering results
% sigma is a fixed covariance for the gaussian mixture.

dim = size(x,2);

% calculate mean of each cluster
mu = zeros(max(idx), size(x,2));
for i=1:max(idx)
    mu(i,:) = mean(x(idx==i,:));
end

% Create GM model
gmm_model = gmdistribution(mu, sigma*ones(1,dim));

% Re assign cluster elements
[gmm_idx] = gmm_model.cluster(x);

end


%% Default parameters
function p = set_default_parameters(varargin)

p = struct;
% Default parameters go here

p.db1_size = 20;
p.db1_minPts = 4;

p.db2_size = 7;
p.db2_minPts = 4;

p.gmm_sigma = 3;

p.ground_truth = [];

p.quiet = 1;

%##########################################################################
% Override default parameters with name/value pairs from the input. 
if numel(varargin) > 0
    p = omex_read_params(p, varargin);
end

end