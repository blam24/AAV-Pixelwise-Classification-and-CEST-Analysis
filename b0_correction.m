function [x_corrected, b0_shift] = b0_correction(x, zspec)
% B0_CORRECTION   Returns B0-corrected frequency offset axis for Z-spectra
%
%   INPUTS:
%       - x, vector of size Mx1 or 1xM containing the prescribed frequency offsets in ppm.
%       - zspec, matrix MxN where each column represents a Z-spectrum.
%
%   OUTPUTS:
%       - x_corrected, matrix MxN where each column represents a B0-corrected frequency offset.
%       - b0_shift, vector Nx1 where each value represents the B0 shift for each of the N Z-spectra.
%
%   USAGE:
%       - [x_corrected, b0_shift] = b0_correction(x, zspec);
%
%   AUTHOR:
%       - Kevin Godines (kevingodines@berkeley.edu)
%
%   DATE:
%       - 2019/10/08

%      M0 |~Water~~~~~~~~~~~~~~~~~~~~~|   |~MT~~~~~~~~~~~~~~~~~~~~~~~~|
%         | Amplitude | FWHM | Center |   | Amplitude | FWHM | Center |
P0 = [ 1    0.8         1.8    0            0.15        40     -1     ];
lb = [ 1    0.02        0.3    -10          0.0         30     -2.5   ];
ub = [ 1    1           10     10           0.5         60     0      ];

% fit-options
options   = [1E-04, 1E-15, 1E-10, 1E-04, 1E-06];
nIter     = 50;
matlab_options = optimset('TolFun',options(4),'TolX',options(3), 'MaxIter',nIter,'Display','off');


n_zspec = size(zspec,2);
x_corrected = zeros(size(zspec));
n_offset = length(x);
b0_shift = zeros(n_zspec,1);
x = reshape(x,[n_offset,1]);

for i = 1:n_zspec % For each Z-spectrum
    %Lorentzian function with CEST pools
    fun1 = @(P,x)P(1) - P(2)*P(3)^2./4./(P(3)^2./4+(x-P(4)).^2)... % Water
        - P(5)*P(6)^2./4./(P(6)^2./4+(x-P(7)).^2);... % MT    
    T = lsqcurvefit(fun1,P0,x,zspec(:,i),lb,ub,matlab_options); %Fitting algorithm
    b0_shift(i) = T(4); %This the B0 offset in ppm obtained from lorentzian fitting.
    x_corrected(:,i) = x - b0_shift(i);
end
