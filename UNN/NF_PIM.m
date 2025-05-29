

clear all

global e0 mu0 c Z0

tic

e0 = 8.8541878128E-12; % electric constant (permittivity of free space) 
mu0 = 1.25663706212E-6; % magnetic constant (permeability of free space)

c = 2.99792458 * 1E8; % % speed of light in vacuum (c = (e0*mu0)^(-1/2) )
Z0 = 376.730313668; % % impedance of free space (Z0 = pi * 119.9169832, Z0 = mu0 * c )

f_0 = 2.6E9;    
f_1 = 2.630E9;
f_2 = 2.680E9;
f_3 = 2 * f_1 - f_2;

% k = 2 * pi * f / c;
% lambda = c ./ f;


% % call the function  that sets parameters of a single antenna element, antenna array and a PIM sources
[ SingleAntennaParameters, AntennaArrayParameters, PIMsourceParameters ] = Parameters();



% % calculate voltage, induced by antenna array elements at PIM sources
% at f_1 and f_2 arrays Voltage_f1,Voltage_f2), and mutual impedances 
% of antenna elements and  PIM sources at  f_3 (array MutualImpedaces_f3 )

% The arrays Voltage_f1, Voltage_f2, MutualImpedaces_f3 have size [N_PIM x N_array], 
% N_PIM - number of PIM sources, 
%     N_array - number of antenna array elements

[ Voltage_f1, Voltage_f2, MutualImpedaces_f3 ] = VoltagesInducedOnPIM( f_1, f_2, f_3,...
    SingleAntennaParameters, AntennaArrayParameters, PIMsourceParameters );


%%---%%---%%---%%---%%-- your code --%%---%%---%%---%%

    % % weight vectors (size [N_array x 1]) 
    % amplitude and phase distribution at f_1 and f_2. 
    % Arrays w_f1 and w_f1 should have size [N_array x 1]. 
    % Normalization of vectors w_f1, w_f1 is performed in function 'PIM3_levels'

    w_f1 = ones(AntennaArrayParameters.N_array, 1); w_f1 = w_f1./sqrt(w_f1'*w_f1);
    w_f2 = ones(AntennaArrayParameters.N_array, 1); w_f2 = w_f2./sqrt(w_f2'*w_f2);

%%---%%---%%---%%---%%-- your code  --%%---%%---%%---%%


% %   
% the total power of one Sc have to be equal to 1
numSc = 600;
AntennaArrayParameters.w1 = bsxfun(@times,w_f1,complex(randn(1,numSc),randn(1,numSc))/sqrt(2)); % weight vector at the frequency f1 (numAnt x numSc)
AntennaArrayParameters.w2 = bsxfun(@times,w_f2,complex(randn(1,numSc),randn(1,numSc))/sqrt(2)); % weight vector at the frequency f2 (numAnt x numSc)

% % calculate voltage at antenna element's loads at f_3
    % PIM3 - has size[N_array x 1], contains normilized voltage values 
[ PIM3_voltages ] = PIM3_levels( f_1, f_2, Voltage_f1, Voltage_f2, MutualImpedaces_f3, ...
    SingleAntennaParameters, AntennaArrayParameters, PIMsourceParameters );


plot(abs(PIM3_voltages(1,:)));


toc

