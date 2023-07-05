% Chamber and Nozzle sizing

% Parameters
chamber_pressure = 180; % PSI; chosen
exit_pressure = 14.7; % PSI; chosen
OF = 3; % Dimensionless; chosen
F = 300; % lbf; chosen
contraction = 4; % dimensionless; chosen
L_star = 1; % m; chosen

% From RPA
T1 = 2524.6864; % Kelvin
rho1 = 1.0352; % kg/m^3
gamma = 1.2291; % dimensionless; ratio of specific heats
R = 368.1; % J/kg-K; gas constant

% Unit Conversions
chamber_pressure = chamber_pressure * 6894.76; % Pa
exit_pressure = exit_pressure * 6894.76; % Pa
F = F * 4.448; % Newtons

%% CALCULATIONS -----------------------------------------------------------

% Find temp at throat (critical temp, Tt)
Tt = (2 * T1) / (gamma + 1); % K

% Find throat velocity (At mach 1, so really the sonic velocity)
vt = (gamma * R * Tt) ^ 0.5; % m/s

% Find exit velocity
v2 = (((2 * gamma) / (gamma - 1)) * R * T1 * (1 - (exit_pressure / chamber_pressure) ^ ((gamma - 1) / gamma))) ^ 0.5; % m/s

% Find mass flow rate
mdot = F / v2; % kg/s

% Find specific volume in chamber/nozzle inlet
V1 = (R * T1) / chamber_pressure; % m^3/kg
rho1_inlet = 1 / V1; % kg/m^3

% Find specific volume at throat
Vt = V1 * ((gamma + 1) / 2) ^ (1 / (gamma - 1)); % m^3/kg
rhot = 1 / Vt; % kg/m^3

% Find specific volume at nozzle exit
press_ratio = (chamber_pressure / exit_pressure) ^ ((gamma - 1) / gamma);
V2 = V1 * (press_ratio) ^ (1 / (gamma - 1)); % m^3/kg
rho2 = 1 / V2; % kg/m^3

% Find throat area
At = (mdot * Vt) / vt; % m^2
Rt = ((At / pi) ^ 0.5) * 1000; % mm

% Find nozzle exit area
Ae = (mdot * V2) / v2; % m^2
Re = ((Ae / pi) ^ 0.5) * 1000; % mm

% Find chamber volume
chamber_volume = L_star * At; % m^3

% Find chamber cross sectional area
chamber_area = contraction * At; % m^2
Rc = ((chamber_area / pi) ^ 0.5) * 1000; % mm

% Find chamber length 
chamber_length = chamber_volume / chamber_area; % m


