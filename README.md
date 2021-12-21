# bGRID_geo
% -------------------------------------------------------------------------
% results = bGRID(catalog_ZMAP, polygon_coord, grid_size, Nmin_b, bin,...
%                   Mc_type)
% Inputs:
%       catalog_ZMAP: catalog with ZMAP format
%       polygon_coord: coordenates of the polygon where the calculations
%                       [min_x max_x; min_y max_y]
%       grid_size: size of the grid [size_x size_y]
%       Nmin_b: minimum number of events to calculate the b-value in each
%       grid cell
%       bin: selection of the bin size for the construction of the
%       frequency-magnitud distribution. A value of 0.1 is recommended.
%       Mc_type: selection of magnitude of completeness calculation
%       methodology. Options are maximum curvature (mc_maxc) and goodness
%       of fit.

% Outputs: 7 or 8 figures and one cell array (results) with the results:
%          Figures:
%           1. Magnitude of completeness vs Goodness of fit (R_GoF) or 
%              Histogram of frequency-magnitude distribution (MaxC). 
%           2. b-value grid plot
%           3. b-value grid plot with seismic events
%           4. Standard error of b-value estimation, sigma_b (Shi & Bolt, 1982)
%           5. Mc_GoF + R_GoF or Mc_MaxC grid plot
%           6. Count above Mc grid plot (Number of events with magnitude
%           greater than Mc per grid)
%           7. Number of events per grid plot (Total number of events per 
%           grid)

%           results cell array {10x2} or {8x2} (per line):
%           1. [1 x n] matrix with magnitudes of the n seismic events inside
%           the polygon with coordinates polygon_coord
%           2. Value of Mc_GoF or Mc_MaxC for all the data inside the polygon.
%           3. R_GoF (only for mc_gof methodology) with bins
%           4. [4 x 1] matrix with global G-R parameters: a, b, sigma_b, R]
%           5. Number of events above Mc
%           6. Expected maximum magnitude, M_max = a/b
%           7. {4 x 1} cell array with grid values of G-R parameters
%           8. [(max_x-min_x)/size_x x (max_y-max_y)/size_y] matrix with
%              number of events above Mc per grid
%           9. [(max_x-min_x)/size_x x (max_y-max_y)/size_y] matrix with
%              the value of Mc per grid
%           10.[(max_x-min_x)/size_x x (max_y-max_y)/size_y] matrix with
%              the value of R_gof per grid (only with mc_gof)

% Authors: 
% Ph. D. Rodrigo Estay Huidobro
%        Universidad Técnica Federico Santa María
%        Santiago, Chile
% Ph. D. Claudia Pavez Orrego
%        Norges Geologiske UndersФkelse
%        Trondheim, Noruega

% Version 1: December/21
%-------------------------------------------------------------------------
