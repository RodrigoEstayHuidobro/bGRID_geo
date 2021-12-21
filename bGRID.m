function results = bGRID(catalog_ZMAP, polygon_coord, grid_size, Nmin_b,...
    bin, Mc_type)
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

%% Variables initiation

long = catalog_ZMAP(:,1);
lat = catalog_ZMAP(:,2);

%% Polygon selection

% Global coordinates / South Norway / Ridge / North Norway / Central / Artic
% Polygon coordinates

min_long = polygon_coord(1); % long min: -10 / 0 / -12 / 7 / 5 / -12
max_long = polygon_coord(3); % long max: 40 / 12 / 25 / 30 / 15 / 25

min_lat = polygon_coord(2); % lat min: 50 / 56 / 71 / 67 / 64 / 71
max_lat = polygon_coord(4); % lat max: 90 / 64 / 79 / 72 / 67 / 79

dlong = grid_size(1);
dlat = grid_size(2);

% MidNor: latmin = 64; latmax = 72
         %longmin = 3.5; longmax = 23.3
         
% Northern: latmin = 68.3; latmax = 72
           %longmin = 19.8; longmax = 34
           
%% Catalog inside polygon

IN_polygon = inpolygon(long, lat, polygon_coord(1,:), polygon_coord(2,:));

catalog_in_polygon = catalog_ZMAP(IN_polygon,:);

long_polygon = catalog_in_polygon(:,1);
lat_polygon = catalog_in_polygon(:,2);
mag_polygon = catalog_in_polygon(:,6);

%% Selection of Mc

switch lower(Mc_type)   
    case 'mc_gof'
        
        [Mc_sector, Rgof_sector, bins_sector, resume_sector] =...
            goodness_of_fit(mag_polygon, bin);
               
        figure
        plot(bins_sector(1:end-1), resume_sector(1:end-1,5))
        ylabel('Rgof')
        xlabel('Magnitude of completeness, Mc')
        grid on
        saveas(gcf,'Rgof_vs_Mc','epsc')
        
    case 'mc_maxc'
        
        [Mc_sector, mag_bin, bins] = Mc_MaxC(mag_polygon, bin);
        
        figure
        hist(mag_bin, bins);
        xlabel('Magnitude')
        ylabel('Frequency')
        
        saveas(gcf,'Histogram','epsc')        
end

%% Main

ind_i = 0;
ind_j = 0;

lat_vector = min_lat:dlat:max_lat-dlat;
long_vector = min_long:dlong:max_long-dlong;

a = zeros(length(lat_vector), length(long_vector));
b = NaN(length(lat_vector), length(long_vector));
sigma_b = zeros(length(lat_vector), length(long_vector));
R = zeros(length(lat_vector), length(long_vector));
Mc_cell = NaN(length(lat_vector), length(long_vector));
Rgof_cell = NaN(length(lat_vector), length(long_vector));
count_above_Mc = zeros(length(lat_vector), length(long_vector));
count_data = zeros(length(lat_vector), length(long_vector));

[a_all, b_all, sigma_b_all, R_all, ~, ~] = GR(mag_polygon,Mc_sector,bin);
numbers_of_events_aboveMc = sum(mag_polygon >= Mc_sector);

Mmax = a_all/b_all;

for i = lat_vector
    
    ind_i = ind_i + 1;
    
    for j = long_vector
        
        ind_j = ind_j + 1;
        
        IN = inpolygon(long_polygon,lat_polygon,[j, j+dlong],[i, i+dlat]);
        
        count_above_Mc(ind_i, ind_j) = sum(mag_polygon(IN) >= Mc_sector);
        count_data(ind_i,ind_j) = sum(IN);
        
        if count_above_Mc(ind_i, ind_j) >= Nmin_b
                        
            try
            [a(ind_i, ind_j), b(ind_i, ind_j), sigma_b(ind_i, ind_j),...
                R(ind_i, ind_j), ~, ~] = GR(mag_polygon(IN),Mc_sector,bin);  
                
            catch ME
                a(ind_i, ind_j) = NaN;
                b(ind_i, ind_j) = NaN;
                sigma_b(ind_i, ind_j) = NaN;
                R(ind_i, ind_j) = NaN;
                
            end
            
            switch lower(Mc_type)
                case 'mc_gof'
                    
                    [Mc_cell(ind_i, ind_j), Rgof_cell(ind_i, ind_j), ~, ~] =...
                        goodness_of_fit(mag_polygon(IN), bin);
                    
                case 'mc_maxc'
                    
                    [Mc_cell(ind_i, ind_j), ~, ~] = Mc_MaxC(mag_polygon(IN), bin);
                    
            end
            
        end
        
    end
    
    ind_j = 0;
    
end

a = flipud(a);
b = flipud(b);
sigma_b = flipud(sigma_b);
R = flipud(R);
count_above_Mc = flipud(count_above_Mc);
Mc_cell = flipud(Mc_cell);
Rgof_cell = flipud(Rgof_cell);

switch lower(Mc_type)
    case 'mc_gof'
        
        results = {'Magnitudes' mag_polygon;...
            'Mc_{GoF} polygon' Mc_sector; 'R_{GoF} polygon' Rgof_sector;...
            'G-R parameters all' [a_all, b_all, sigma_b_all, R_all];...
            'Number of events above Mc' numbers_of_events_aboveMc;...
            'Expected M_{max}' Mmax;...
            'G-R parameters grid' {a; b; sigma_b; R;};...
            'Count above Mc' count_above_Mc;...
            'Mc grid' Mc_cell; 'R_{GoF} grid' Rgof_cell};
        
    case 'mc_maxc'
        
        results = {'Magnitudes' mag_polygon;...
            'Mc_{MaxC} polygon' Mc_sector;...
            'G-R parameters all' [a_all, b_all, sigma_b_all, R_all];...
            'Number of events above Mc' numbers_of_events_aboveMc;...
            'Expected M_{max}' Mmax;...
            'G-R parameters grid' {a; b; sigma_b; R;};...
            'Count above Mc' count_above_Mc;...
            'Mc grid' Mc_cell};
        
end

%% Figures

figure
imagesc(long_vector+dlong/2,lat_vector+dlat/2,flipud(b))
set(gca,'YDir','normal')
caxis([0 max(max(b))])
colorbar
colormap([1 1 1; parula(256)])
xlabel('Longitude [°]')
ylabel('Latitude [°]')
title('b-value')
saveas(gcf,'b','epsc')

figure
imagesc(long_vector+dlong/2,lat_vector+dlat/2,flipud(b))
set(gca,'YDir','normal')
caxis([0 max(max(b))])
colorbar
colormap([1 1 1; parula(256)])
xlabel('Longitude [°]')
ylabel('Latitude [°]')
title('b-value + seismic events')
hold on
plot(long_polygon, lat_polygon, 'r*')
saveas(gcf,'b + seismic events','epsc')

figure
imagesc(long_vector+dlong/2,lat_vector+dlat/2,flipud(sigma_b))
set(gca,'YDir','normal')
caxis([0 max(max(sigma_b))])
colorbar
colormap([1 1 1; parula(256)])
xlabel('Longitude [°]')
ylabel('Latitude [°]')
title('Estimation error \sigma_b')
saveas(gcf,'sigma_b','epsc')

switch lower(Mc_type)
    case 'mc_gof'
        
        figure
        imagesc(long_vector+dlong/2,lat_vector+dlat/2,flipud(Mc_cell))
        set(gca,'YDir','normal')
        caxis([0 max(max(Mc_cell))])
        colorbar
        colormap([1 1 1; parula(256)])
        xlabel('Longitude [°]')
        ylabel('Latitude [°]')
        title('Mc GoF')
        saveas(gcf,'Mc_GoF','epsc')
        
        figure
        imagesc(long_vector+dlong/2,lat_vector+dlat/2,flipud(Rgof_cell))
        set(gca,'YDir','normal')
        caxis([0 max(max(Rgof_cell))])
        colorbar
        colormap([1 1 1; parula(256)])
        xlabel('Longitude [°]')
        ylabel('Latitude [°]')
        title('Rgof')
        saveas(gcf,'Rgof','epsc')
        
    case 'mc_maxc'
        
        figure
        imagesc(long_vector+dlong/2,lat_vector+dlat/2,flipud(Mc_cell))
        set(gca,'YDir','normal')
        caxis([0 max(max(Mc_cell))])
        colorbar
        colormap([1 1 1; parula(256)])
        xlabel('Longitude [°]')
        ylabel('Latitude [°]')
        title('Mc MaxC')
        saveas(gcf,'Mc_MaxC','epsc')
        
end

figure
imagesc(long_vector+dlong/2,lat_vector+dlat/2,flipud(count_above_Mc))
set(gca,'YDir','normal')
caxis([0 max(max(count_above_Mc))*.8])
colorbar
colormap([1 1 1; parula(256)])
xlabel('Longitude [°]')
ylabel('Latitude [°]')
title('Count above Mc')
saveas(gcf,'Count_above_Mc','epsc')

figure
imagesc(long_vector+dlong/2,lat_vector+dlat/2,count_data)
set(gca,'YDir','normal')
caxis([0 max(max(count_data))])
colorbar
xlabel('Longitude [°]')
ylabel('Latitude [°]')
title('Number of events per grid')
saveas(gcf,'Number_of_events_per_grid','epsc')

end

function [Mc, mag_bin, bins] = Mc_MaxC(Mw, bin)

mag_bin = mag2bin(Mw, bin);

bins = min(mag_bin):bin:max(mag_bin);

N = hist(mag_bin, bins);

Mc = bins(find(N == max(N),1));

end

function [Mc, Rgof_final, bins, resume] = goodness_of_fit(data, bin)

M = mag2bin(data, bin);

bins = min(M):bin:max(M);

a = zeros(length(bins),1);
b = zeros(length(bins),1);
sigma_b = zeros(length(bins),1);
R_GR = zeros(length(bins),1);
Rgof = zeros(length(bins),1);
N = zeros(length(bins));

frec = hist(M,bins)';

B = zeros(length(frec),1);

for k = 1:length(frec)
    
    B(k,1) = sum(frec(k:end));
    
end

for i = 1:length(bins)-1
    
    Mci = round(bins(i)*10)/10;
    
    [a(i), b(i), sigma_b(i), R_GR(i), N(i:end,i)] = GR(M, Mci, bin);
 
    Rgof(i) = 100 - sum(abs(B(i:end)-N(i:end,i)))/sum(B(i:end))*100;
    
end

Rgof_final = Rgof(find(Rgof >= 90, 1));

Mc = bins(find(Rgof >= 90, 1));

if isempty(Mc)
    
    Rgof_final = max(Rgof);
    
    Mc = bins(Rgof == max(Rgof));

end

resume = [a, b, sigma_b, R_GR, Rgof];

end

function mag_bin = mag2bin(mag,bin)

% Función que pasa de magnitudes decimales (MAG) a magnitudes en bins (BIN)
% aproximadas en nPlaces decimales.

nPlaces = length(num2str(bin)) - 2;

fneg = find(mag < 0);
fpos = find(mag > 0);

mag(fneg) = -roundTo(-mag(fneg)-bin/100,nPlaces);
mag(fpos) = roundTo(mag(fpos),nPlaces);

mag_bin = mag;

end

function [a, b, sigma_b, R, Ncalc, AIC] = GR(data,Mc,bin)

% Cálculo de b y R (Gutenberg - Richter).
% Función que calcula los parámetros de Gutenberg-Richter para un vector de
% magnitudes DATA y una magnitud de completitud Mc. BIN corresponde al
% tamaño de los bins que se quiere utilizar para agrupar las magnitudes. Se
% recomienda BIN = 0.1 (Marzocchi & Sandri, 2003).
% Se asume que la magnitud sigue una densidad de probabilidad exponencial
% con media 1/beta [f = beta*exp(-beta*(M-Mc))]

%% Inicialización de variables y condiciones iniciales

Mmax = max(data);

KM = floor(10*(Mmax-Mc+bin/10));

if KM == 0
    
    a = NaN;
    b = NaN;
    sigma_b = NaN;
    R = NaN;
    Ncalc = NaN;
    AIC = NaN;
    
else 
    %% Cálculo de a y b
    
    K = KM:-1:0;
      
    bins = Mc:bin:Mmax;
    
    frec = hist(data,bins)';
    frec(1) = frec(1) - length(find(data < Mc));
        
    Nacum = zeros(length(frec),1);
    
    for k = 1:length(frec)
        
        Nacum(k,1) = sum(frec(k:end));
        
    end
    
    NC = max(Nacum);
    
    NK = flipud(frec);
    NK(:,2) = cumsum(NK(:,1));
    
    aux = 0;
    SM = zeros(length(NK),1);
    
    for k = 1:length(NK)
        
        SM(k,1) = (bin*K(k)+bin/2)*NK(k)+aux;
        aux = SM(k,1);
        
    end
    
    SMmax = max(SM);
    
    beta = NC/SMmax;
    
    b = beta/log(10);
    
    Mz = Mc - bin/2;
    %M1 = log(NC)/beta+Mz;
    a = (log(beta*NC) + beta*Mz)/log(10);
    
    %% Cálculo de sigma_b (Shi & Bolt, 1982)
    M_mean = 1/beta + Mc - bin/2;
    dif_cuad = (bins - M_mean).^2;
    sigma_b = log(10)*b^2*sqrt(sum(dif_cuad)/(NC*(NC-1)));
    
    %% Cálculo de R
    Ncalc = 10.^(a - b.*(bins-bin/2))./beta;
    Ncalc = Ncalc';
    R = 1 - sum(abs(Ncalc - Nacum))/sum(Ncalc);
    
    %     figure
    %     plot(bins,log10(Nacum),'*',bins,log10(Ncalc),'-r')
    
    %% AIC
    
    AIC = -2*NC*(log(beta)-1)+2;
    
end

end