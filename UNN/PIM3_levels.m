function [ PIM3 ] = PIM3_levels( f_1, f_2, Voltage_f1, Voltage_f2, MutualImpedaces_f3, ...
    SingleAntennaParameters, AntennaArrayParameters, PIMsourceParameters )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global c 
%[ E_x3, E_y3, E_z3, H_x3, H_y3, H_z3 ] = PIM_response( x,y,z, f_1, f_2, x_PIM, y_PIM, z_PIM, E_xyz1, E_xyz2, PIM_parameters )


% % radiated power at f_1, f_2
Power_f1 = 1E-3 * 10^(AntennaArrayParameters.Power_f1 / 10);
Power_f2 = 1E-3 * 10^(AntennaArrayParameters.Power_f1 / 10);

% % weight vectors
w_f1 = AntennaArrayParameters.w1; % weight vector at the frequency f1
w_f2 = AntennaArrayParameters.w2; % weight vector at the frequency f2
% % we normalize the weight vectors
numSc1 = size(w_f1,2);
numSc2 = size(w_f2,2);
numPim = size(Voltage_f1,1);
w_f1 = sqrt(Power_f1) * w_f1 ./sqrt(numSc1); % numAnt x numSc
w_f2 = sqrt(Power_f2) * w_f2 ./sqrt(numSc2);

% % Calculate voltage at the PIM source taking into account amplitude an d phase distribution at the Antenna Array
V_f1 = Voltage_f1 * w_f1; %numPim x numSc
V_f2 = Voltage_f2 * w_f2;

% % % calculate current at f_3
a1 = PIMsourceParameters.PolynomNonlinearityCoeff_a1;
a3 = PIMsourceParameters.PolynomNonlinearityCoeff_a3;

Z_L = 1 ./ a1;

L_PIM = PIMsourceParameters.L_PIM; %m, dipole antenna length
a_PIM = PIMsourceParameters.a_PIM;


Z_inp1 = 20 * pi^2 * (f_1 .* L_PIM / c).^2 - 1i * 120 * ( log(L_PIM./(2*a_PIM)) - 1 ) ./ tan(pi * f_1 .* L_PIM / c);
Z_inp2 = 20 * pi^2 * (f_2 .* L_PIM / c).^2 - 1i * 120 * ( log(L_PIM./(2*a_PIM)) - 1 ) ./ tan(pi * f_2 .* L_PIM / c);

V_inp1 = bsxfun(@times,V_f1, Z_L ./ (Z_inp1 + Z_L)); %numPim x numSc
V_inp2 = bsxfun(@times, V_f2, Z_L ./ (Z_inp2 + Z_L));

% I_f3 = (3/4) * a3 .* V_inp1.^2 .* conj(V_inp2); % single carrier

% WB signal processing
I_f3 = zeros(numPim,2*numSc1+numSc2-2);
for cnPim = 1:numPim
    A = V_inp1(cnPim,:).'*V_inp1(cnPim,:);
    A = A(:,end:-1:1);
    b = NaN(2*numSc1-1,1);
    for k = -numSc1+1:numSc1-1
        b(k+numSc1) = sum(diag(A,k));
    end
    I_f3(cnPim,:) = (3/4)*a3(cnPim)*conv(b,V_inp2(cnPim,:)');
    I_f3(cnPim,:) = I_f3(cnPim,end:-1:1);
end

% % Calculate voltage  at the Antenna Array ellements using  matrix of mutual impedances 
VoltagesOnArray_f3 = MutualImpedaces_f3.' * I_f3;


% % calculate voltage at antenna matched loads
VoltagesOnAntennaLoad = VoltagesOnArray_f3 / 2;

% % normalized load
NormalizedVoltagesOnAntennaLoad = VoltagesOnAntennaLoad / sqrt(real(SingleAntennaParameters.Z0));

% PIM3 = 10 * log10(abs(VoltagesOnArray_f3) / 1E-3);
% PIM3 = 10 * log10(abs(VoltagesOnAntennaLoad) / 1E-3);

PIM3 = NormalizedVoltagesOnAntennaLoad;

end

