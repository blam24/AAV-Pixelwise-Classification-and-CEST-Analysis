function F = zspec_lorentzian_fit_hydrogel(x_corrected, zspec, method)
% ZSPEC_LORENTZIAN_FIT   Fits z-spectra to a sum of Lorentzians using a variety of methods. 
%
%   INPUTS:
%       - x_corrected, matrix MxN where each column N represents a B0-corrected frequency offset.
%       - zspec, matrix MxN where each column represents a Z-spectrum.
%       - method, string or char array for fitting method.
%                '2stepLD' for 2-step lorentzian difference
%                '1step' for single-step fitting.
%
%   OUTPUTS:
%       - F, struct containing:
%           - F.zspec_fitted
%           - F.x_interp
%           - F.fit_parameters
%           - F.water
%           - F.mt
%           - F.amine
%           - F.noe
%           - F.aav
%           - F.amide
%           - F.zspec
%           - F.lorentzian_difference
%           - F.x_corrected
%           - F.nmse.full_zspec
%           - F.nmse.aav
%           - F.nmse.amide
%           - F.contrasts.water
%           - F.contrasts.amine
%           - F.contrasts.noe
%           - F.contrasts.aav
%           - F.contrasts.amide
%
%   USAGE:
%       - method = '2stepLD'; % alternatively: method = '1step';
%         F = zspec_lorentzian_fit(x_corrected, zspec, method);
%
%   AUTHOR:
%       - Bonnie Lam (bonnie.lam@berkeley.edu)
%       - Kevin Godines (kevingodines@berkeley.edu)
%
%   DATE:
%       - 2022/10/25

%% Setting up variables
n_seg = size(zspec,2);
n_interp = 1000;
cutoffs = [-3, -1.4, 1.4, 3]; % exclude points between -3ppm and -1.4 ppm, and between 1.4ppm and 3ppm (CEST metabolites region)

d1 = zeros(1,n_seg);
d2 = zeros(1,n_seg);
d3 = zeros(1,n_seg);
d4 = zeros(1,n_seg);
ix1 = zeros(1,n_seg);
ix2 = zeros(1,n_seg);
ix3 = zeros(1,n_seg);
ix4 = zeros(1,n_seg);

for i = 1:n_seg
    [d1(i), ix1(i)] = min(abs(x_corrected(1:ceil(size(x_corrected,1)/2),i)-(cutoffs(1))));
    [d2(i), ix2(i)] = min(abs(x_corrected(1:ceil(size(x_corrected,1)/2),i)-(cutoffs(2))));
    [d3(i), ix3(i)] = min(abs(x_corrected(ceil(size(x_corrected,1)/2):end,i)-(cutoffs(3))));
    ix3(i) = ix3(i) + ceil(size(x_corrected,1)/2) -1;
    [d4(i), ix4(i)] = min(abs(x_corrected(ceil(size(x_corrected,1)/2):end,i)-(cutoffs(4))));
    ix4(i) = ix4(i) + ceil(size(x_corrected,1)/2) -1;
    x_cropped{i} = [x_corrected(1:ix1(i),i);x_corrected(ix2(i):ix3(i),i); x_corrected(ix4(i):end,i)]; % Setting cropped frequency axis
    x_interp(:,i) = linspace(x_corrected(1,i),x_corrected(end,i),n_interp); % Defining frequency offset axis for smooth interpolation
    zspec_cropped{i} = zspec([1:ix1(i),ix2(i):ix3(i),ix4(i):end],i); %cropped z-spectra
end

clear d1 d2 d3 d4 ix1 ix2 ix3 ix4

% fit-options
options   = [1E-04, 1E-15, 1E-10, 1E-04, 1E-06];
nIter     = 50;
matlab_options = optimset('TolFun',options(4),'TolX',options(3), 'MaxIter',nIter,'Display','off');

if not(exist('method','var'))
    method = '1step';
    warning('Method not found. Defaulting to single step fitting');
end

if strcmp(method,'2stepLD')
    %% Lorentzian fitting: LD 2-pool model (background)
    % parameter order for each pool is: amplitude, FWHM, peak center
    %      [M0]  |Water~~~~~~~~~|   |MT~~~~~~~~~~~|
    P0 = [ 1    0.8    0.8     0     0.02  40   -1]; % starting point
    lb = [ 1    0.02   0.2     0     0.0   30   -2.5]; % lower bound
    ub = [ 1    1.0    5.0     0     0.15  60    0]; % upper bound
    
    fit_parameters_2pool = cell(n_seg,1);
    background = zeros(size(x_corrected));
    water = zeros(size(x_interp));
    mt = zeros(size(x_interp));
    
    for i = 1:n_seg %for each segment
        fun1 = @(P,x)P(1) - P(2)*P(3)^2./4./(P(3)^2./4+(x-P(4)).^2)... % Water
            - P(5)*P(6)^2./4./(P(6)^2./4+(x-P(7)).^2); % MT
        
        Q = lsqcurvefit(fun1,P0,x_cropped{i},zspec_cropped{i},lb,ub,matlab_options); %Fitting algorithm using cropped freq axis
        fit_parameters_2pool{i} = Q;
        background(:,i) =  Q(2)*Q(3)^2./4./(Q(3)^2./4+(x_corrected(:,i)-Q(4)).^2) + Q(5)*Q(6)^2./4./(Q(6)^2./4+(x_corrected(:,i)-Q(7)).^2);
        water(:,i) = Q(2)*Q(3)^2./4./(Q(3)^2./4+(x_interp(:,i)-Q(4)).^2);
        mt(:,i) = Q(5)*Q(6)^2./4./(Q(6)^2./4+(x_interp(:,i)-Q(7)).^2);
    end
    
    %% Generate Lorentzian Difference Contrast
    lorentzian_difference = zeros(size(x_corrected));
    for i = 1:n_seg % for each segment
        lorentzian_difference(:,i) = 1-(zspec(:,i) + background(:,i));
    end
    
    %% Lorentzian fitting: LD 4-pool model (metabolites)
    % parameter order for each pool is: amplitude, FWHM, peak center
    %     % |NOE~~~~~~~~~~~~~|  |Amine~~~~~~~~~~~~|  |Amide~~~~~~~~~~~~|  |AAV~~~~~~~~~~~~|
    P0 = [ 	0.05   3    -2.75    0.025  0.7   1.8	  0.025  0.7   3.6     0.03   0.7  0.7];  % starting point
    lb = [ 	0.0    1    -4.5     0.0    0.4   1.6     0.0    0.4   3.55    0.0    0.4  0.55];  % lower bound
    ub = [ 	0.25   5    -1.5     0.2    1.3   2.0	  0.1    1.3   3.65    0.2    1.3  0.9];  % upper bound
    
    fit_parameters_4pool = cell(1,n_seg);
    fit_parameters = cell(1,n_seg);
    noe = zeros(length(x_interp),n_seg);
    amine = zeros(length(x_interp),n_seg);
    amide = zeros(length(x_interp),n_seg);
    aav = zeros(length(x_interp),n_seg);
    zspec_fitted = zeros(length(x_interp),n_seg);

    noe_nmse = zeros(size(x_corrected,1),n_seg);
    amine_nmse = zeros(size(x_corrected,1),n_seg);
    amide_nmse = zeros(size(x_corrected,1),n_seg);
    aav_nmse = zeros(size(x_corrected,1),n_seg);
    
    for i = 1:n_seg
        % Lorentzian function with CEST pools
        fun2 = @(P,x)P(1)*P(2)^2./4./(P(2)^2./4+(x-P(3)).^2)... % NOE
            + P(4)*P(5)^2./4./(P(5)^2./4+(x-P(6)).^2)... % Amine
            + P(7)*P(8)^2./4./(P(8)^2./4+(x-P(9)).^2)... % Amide
            + P(10)*P(11)^2./4./(P(11)^2./4+(x-P(12)).^2);   % AAV
        
        P = lsqcurvefit(fun2,P0,x_corrected(:,i),lorentzian_difference(:,i),lb,ub,matlab_options); % fitting algorithm
        fit_parameters_4pool{i} = P;
        fit_parameters{i} = [fit_parameters_2pool{i}, fit_parameters_4pool{i}];
        noe(:,i) = P(1)*P(2)^2./4./(P(2)^2./4+(x_interp(:,i)-P(3)).^2);
        amine(:,i) = P(4)*P(5)^2./4./(P(5)^2./4+(x_interp(:,i)-P(6)).^2);
        amide(:,i) = P(7)*P(8)^2./4./(P(8)^2./4+(x_interp(:,i)-P(9)).^2);
        aav(:,i) = P(10)*P(11)^2./4./(P(11)^2./4+(x_interp(:,i)-P(12)).^2);
        zspec_fitted(:,i) = 1 - water(:,i) - mt(:,i) - noe(:,i) - amine(:,i) - amide(:,i) - aav(:,i);
        
        noe_nmse(:,i) = P(1)*P(2)^2./4./(P(2)^2./4+(x_corrected(:,i)-P(3)).^2);
        amine_nmse(:,i) = P(4)*P(5)^2./4./(P(5)^2./4+(x_corrected(:,i)-P(6)).^2);
        amide_nmse(:,i) = P(7)*P(8)^2./4./(P(8)^2./4+(x_corrected(:,i)-P(9)).^2);
        aav_nmse(:,i) = P(10)*P(11)^2./4./(P(11)^2./4+(x_corrected(:,i)-P(12)).^2);
    end
    
elseif strcmp(method,'1step')
    %% Lorentzian fitting: single step 6-pool model
    % parameter order for each pool is: amplitude, FWHM, peak center
    %     |M0| |Water~~~~~~~~~|  |MT~~~~~~~~~~~~~|  |NOE~~~~~~~~~~~~~~|  |Amine~~~~~~~~~~~|  |Amide~~~~~~~~~~~|  |AAV~~~~~~~~~~~~~|
    P0 = [ 1    0.8   1.8   0     0.15  40   -1      0.05   3    -2.75    0.05   0.5   2.0	  0.05   4.5   3.5    0.05   0.5   0.6]; % starting point
    lb = [ 1    0.02  0.3   0     0.0   30   -2.5    0.0    1    -4.5     0.0    0.2   1.6    0.0    0.4   3.2    0.0    0.1   0.55]; % lower bound
    ub = [ 1    1     10    0     0.5   60    0      0.25   5    -1.5     0.3    1     2.4	  0.2    6     3.8    0.2    1.2   0.75]; % upper bound
    
    fit_parameters = cell(1,n_seg);
    water = zeros(size(x_interp));
    mt = zeros(size(x_interp));
    noe = zeros(length(x_interp),n_seg);
    amine = zeros(length(x_interp),n_seg);
    amide = zeros(length(x_interp),n_seg);
    aav = zeros(length(x_interp),n_seg);
    zspec_fitted = zeros(length(x_interp),n_seg);
    
    lorentzian_difference = zeros(size(x_corrected));
    background = zeros(size(x_corrected));
    noe_nmse = zeros(length(x_corrected),n_seg);
    amine_nmse = zeros(length(x_corrected),n_seg);
    amide_nmse = zeros(length(x_corrected),n_seg);
    aav_nmse = zeros(length(x_corrected),n_seg);
    
    for i = 1:n_seg
        % Lorentzian function with CEST pools
        fun3 = @(P,x)P(1) - P(2)*P(3)^2./4./(P(3)^2./4+(x-P(4)).^2)... % Water
            - P(5)*P(6)^2./4./(P(6)^2./4+(x-P(7)).^2)... % MT
            - P(8)*P(9)^2./4./(P(9)^2./4+(x-P(10)).^2)... % NOE
            - P(11)*P(12)^2./4./(P(12)^2./4+(x-P(13)).^2)... % Amine
            - P(14)*P(15)^2./4./(P(15)^2./4+(x-P(16)).^2)... % Amide
            - P(17)*P(18)^2./4./(P(18)^2./4+(x-P(19)).^2);  % AAV
        
        P = lsqcurvefit(fun3,P0,x_corrected(:,i),zspec(:,i),lb,ub,matlab_options); %Fitting algorithm
        fit_parameters{i} = P;
        water(:,i) = P(2)*P(3)^2./4./(P(3)^2./4+(x_interp(:,i)-P(4)).^2);
        mt(:,i) = P(5)*P(6)^2./4./(P(6)^2./4+(x_interp(:,i)-P(7)).^2);
        noe(:,i) = P(8)*P(9)^2./4./(P(9)^2./4+(x_interp(:,i)-P(10)).^2);
        amine(:,i) = P(11)*P(12)^2./4./(P(12)^2./4+(x_interp(:,i)-P(13)).^2);
        amide(:,i) = P(14)*P(15)^2./4./(P(15)^2./4+(x_interp(:,i)-P(16)).^2);
        aav(:,i) = P(17)*P(18)^2./4./(P(18)^2./4+(x_interp(:,i)-P(19)).^2);
        zspec_fitted(:,i) = 1 - water(:,i) - mt(:,i) - noe(:,i) - amine(:,i) - amide(:,i) - aav(:,i);
        
        background(:,i) =  P(2)*P(3)^2./4./(P(3)^2./4+(x_corrected(:,i)-P(4)).^2) + P(5)*P(6)^2./4./(P(6)^2./4+(x_corrected(:,i)-P(7)).^2);
        lorentzian_difference(:,i) = 1-(zspec(:,i) + background(:,i));
        
        noe_nmse(:,i) = P(8)*P(9)^2./4./(P(9)^2./4+(x_corrected(:,i)-P(10)).^2);
        amine_nmse(:,i) = P(11)*P(12)^2./4./(P(12)^2./4+(x_corrected(:,i)-P(13)).^2);
        amide_nmse(:,i) = P(14)*P(15)^2./4./(P(15)^2./4+(x_corrected(:,i)-P(16)).^2);
        aav_nmse(:,i) = P(17)*P(18)^2./4./(P(18)^2./4+(x_corrected(:,i)-P(19)).^2);
    end
    
else
    warning('Fitting failed. Please choose a valid method. Returning empty results.');
    F = struct([]);
    return
end

%% NMSE calculation
fit_all = zeros(size(aav_nmse));
ref_raw = zeros(size(zspec));
fit_nmse_all = zeros(1,n_seg);

for i = 1:n_seg
    fit_all(:,i) = 1 - (background(:,i) + noe_nmse(:,i) + amine_nmse(:,i) + amide_nmse(:,i) + aav_nmse(:,i));
    ref_raw(:,i) = zspec(:,i);
    fit_nmse_all(i) = goodnessOfFit(fit_all(:,i),ref_raw(:,i),'NMSE');
end

%% Storing variables in struct
F = struct();
F.zspec_fitted = zspec_fitted;
F.x_interp = x_interp;
F.fit_parameters = fit_parameters;

F.water = water;
F.mt = mt;
F.noe = noe;
F.amine = amine;
F.amide = amide;
F.aav = aav;
F.zspec = zspec;

if exist('lorentzian_difference','var')
    F.lorentzian_difference = lorentzian_difference;
end

F.x_corrected = x_corrected;

F.nmse.full_zspec = fit_nmse_all;

for i = 1:n_seg
    F.contrasts.water(i) = 100*fit_parameters{i}(2);
    F.contrasts.mt(i) = 100*fit_parameters{i}(5);
    F.contrasts.noe(i) = 100*fit_parameters{i}(8);
    F.contrasts.amine(i) = 100*fit_parameters{i}(11);
    F.contrasts.amide(i) = 100*fit_parameters{i}(14);
    F.contrasts.aav(i) = 100*fit_parameters{i}(17);
end
end