function [ SingleAntennaParameters, AntennaArrayParameters, PIMsourceParameters ] = Parameters( )
% In this function, we set the parameters of the antenna, antenna array and PIM sources

    % SingleAntennaParameters - parameters of a single antenna
    % AntennaArrayParameters - parameters of an antenna array
    % PIMsourceParameters - parameters of PIM sources
    

global e0 mu0 c Z0

c = physconst('LightSpeed');

f = 2.6E9; % GHz 'reference' frequency
lambda = c / f; % wavelength in free spaace at f = f0


%% Parameters of a single antenna

% Choose the type of the single antenna
AntennaType = 'PrintedDipole'; % Antenna type # 1 (Printed dipole)
% AntennaType = 'CylindricalDipole'; % Antenna type # 2 (Cylindrical dipole)

% Here, we set the parameters for each type of antennas
switch AntennaType
    case 'PrintedDipole' % % % Antenna type # 1 (Printed dipole)
        SingleAntennaParameters.AntennaType = 'PrintedDipole_CST'; % % antenna type
        SingleAntennaParameters.CalculationTechnique = 'LinearDipole'; % EM field calculation technique: 'PointDipole' 'LinearDipole' 'Green' 'Multipole' 'Combined_Green+Multipole' 
        SingleAntennaParameters.height = 0.02 * 1E-3; % m (antenna thickness, along the x axis)
        SingleAntennaParameters.width = 4.829 * 1E-3; % m (antenna width, along the y axis)
        SingleAntennaParameters.length = 49.058 * 1E-3; % m (antenna length, along the z axis)
        SingleAntennaParameters.gap = 0.450 * 1E-3; % m (antenna gap, along the z axis)

        SingleAntennaParameters.R0 = 0.5 * lambda; % lower boundary of applicability of the Multipole expansion  
        SingleAntennaParameters.lmax = 5; % maximum multipole order in the Multipole expansion

        SingleAntennaParameters.I0 = 1.1026079E-01 + 1i * 6.9481505E-04; % the input current of an antenna at f = 2.6 GHz (from CST)
        SingleAntennaParameters.Z0 = 8.2245502E+01 + 1i * ( -1.0725507E+00); % the input impedance of an antenna at f = 2.6 GHz (from CST)

        % % loading the file with the electric-current density distribution (from CST)
        load('Current density\Current density - CST_SinglePatch_v.2.mat'); 
        SingleAntennaParameters.Current_density = Current_density; % current density distribution on the antenna

        % calculating the multipole coefficients for the given current density distribution
        [ AE_lm, AM_lm ] = Multipole_Coefficients( f, SingleAntennaParameters); 
        SingleAntennaParameters.AE_lm = AE_lm;
        SingleAntennaParameters.AM_lm = AM_lm;


    case 'CylindricalDipole' % % % % Antenna type # 2 (Cylindrical dipole)
        SingleAntennaParameters.AntennaType = 'CylindricalDipole_CST'; % % antenna type
        SingleAntennaParameters.CalculationTechnique = 'Green'; % EM field calculation technique: 'PointDipole' 'LinearDipole' 'Green' 'Multipole' 'Combined_Green+Multipole' 

        SingleAntennaParameters.width = 5 * 1E-3; % m (antenna diameter)
        SingleAntennaParameters.length = 50 * 1E-3; % m (antenna length, along the z axis)
        SingleAntennaParameters.gap = 0.5 * 1E-3; % m (antenna gap, along the z axis)

        SingleAntennaParameters.R0 = 0.5 * lambda; % lower boundary of applicability of the Multipole expansion  
        SingleAntennaParameters.lmax = 5; % maximum multipole order in the Multipole expansion

        SingleAntennaParameters.I0 = 1.2881514E-01 + 1i * 3.3951765E-02; % the input current of an antenna at f = 2.6 GHz (from CST)
        SingleAntennaParameters.Z0 = 5.2608507E+01 + 1i * (-2.7146733E+01); % the input impedance of an antenna at f = 2.6 GHz (from CST)

        % % loading the file with the electric-current density distribution (from CST)
        load('Current density\Current density - CST_CylindricalDipole_v.1.mat'); 
        SingleAntennaParameters.Current_density = Current_density; % current density distribution on the antenna

        % calculating the multipole coefficients for the given current density distribution
        [ AE_lm, AM_lm ] = Multipole_Coefficients( f, SingleAntennaParameters); 
        SingleAntennaParameters.AE_lm = AE_lm;
        SingleAntennaParameters.AM_lm = AM_lm;

end


%% Parameters of the antenna array

Ny_array = 4; % number of antennas along the y axis
Nz_array = 4; % number of antennas along the z axis
N_array = Ny_array * Nz_array; % total number of antennas in the array

Ly = 0.5 * lambda; % array period along the y axis
Lz = 0.5 * lambda; %  array period along the z axis

% we calculate two-dimensional  of the antenns coordinates of the elements of the antenna array
[Y_array, Z_array] = meshgrid( ([1 : Ny_array] - (Ny_array + 1)/2) * Ly, ([1 : Nz_array] - (Nz_array + 1)/2) * Lz);

% % XYZ_array has size[N_array x 3] and includes x-, y- è z-coordinates of antenna array elements
XYZ_array = zeros(N_array, 3);
XYZ_array(:, 2) = reshape(Y_array,[N_array,1]);
XYZ_array(:, 3) = reshape(Z_array,[N_array,1]);

% % weight vectors of the antenna array
w1 = ones(N_array,1);
w2 = ones(N_array,1);

AntennaArrayParameters.N_array = N_array; % total number of antennas in the antenna array
AntennaArrayParameters.XYZ_array = XYZ_array; % coordinetes of antennas in array
AntennaArrayParameters.w1 = w1 / sqrt(w1' * w1); % weight vector at the frequency f1
AntennaArrayParameters.w2 = w2 / sqrt(w2' * w2); % weight vector at the frequency f2
AntennaArrayParameters.Power_f1 = 40; % dBm, total radiated power at the frequency f1
AntennaArrayParameters.Power_f2 = 40; % dBm, total radiated power at the frequency f2


%% Parameters of the PIM sources

%%--%%--%%--%%--%%--%%--%% 
% Code between lines %%--%%, corresponds to sertain set of parameters of PIMs. 
% One can comment it and copy to get another set of parameters

    N_PIM = 3; % total number of the PIM sources

    % setting the coordinates of PIM sources
    % lines correspondes to different PIMs, columns- x-,y- and z- coordinates
    XYZ_PIM = [1 0 -1; 
               1 0 0; 
               1 0 1] * lambda ; % %coordinates of PIM sources

    % setting the orientation of PIM sources 
    % lines correspondes to different PIMs, columns- x-,y- and z- projection of orientation vector
    % (PIM_orientation())
    PIM_orientation = [1 0 0; 
                       0 1 0; 
                       0 0 1]; % vectors determing the orientation of PIM sources   

    % setting the length and radius of the PIM-source antenna               
    L_PIM = [10; 15; 20] * 1E-3; %m, PIM antenna length
    a_PIM = [1; 1; 1] * 1E-3; %m, PIM antenna radius

    % Parameters of the Polynomial model of a nonlinear Current-Voltage function for the antenna load
        Vt = 26 * 1E-3; % thermal voltage in the diode I-V characteristic
        Is = 1E-6; % the saturation current in the diode I-V characteristic
        I0 = [1; 1; 1] * Is;
        V0 = [1; 1; 1;] * Vt;
        a1 = 2 * I0 ./ V0; %  
        a3 = I0 ./ (3 * V0.^3);
%%--%%--%%--%%--%%--%%--%%              

PIMsourceParameters.N_PIM = N_PIM;% total number of the PIM sources
PIMsourceParameters.XYZ_PIM = XYZ_PIM; % coordinates of PIM sources
PIMsourceParameters.PIM_orientation = PIM_orientation; % vectors determing the orientation of PIM sources
PIMsourceParameters.L_PIM = L_PIM; %m, PIM antenna length
PIMsourceParameters.a_PIM = a_PIM; %m, PIM antenna radius
PIMsourceParameters.PolynomNonlinearityCoeff_a1 = a1;
PIMsourceParameters.PolynomNonlinearityCoeff_a3 = a3;



end

