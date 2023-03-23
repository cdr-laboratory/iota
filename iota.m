%
% #########################################################################
% ### iota -- a stochastic model of the marine particle factory ###########
% #########################################################################
%
% -------------------------------------------------------------------------
% NOTE: based on Fakhraee et al. [2023]
% NOTE: preprint doi:10.1038/s41561-020-006600-6
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- CURATE VARIABLES AND PATHS ------------------------------------------
% -------------------------------------------------------------------------
clear all
tic
% -------------------------------------------------------------------------
addpath src/;           % make sure this points to src directory
addpath input/;    
load calibration_final
% -------------------------------------------------------------------------
global F_NPP F_mineral velocity_total z mineral_mass
global Kz_water Adv_TH Adv_TH1 Rate_org_interp TAinit Adv_TH2 
global kO2 Iron O2init  DICinit kFe R1_carb R_diss_mineral_1
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- SET INPUT PARAMETERS ------------------------------------------------
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- material fluxes and silicate mineral options ------------------------
% -------------------------------------------------------------------------
F_mineral        = 0;          % silicate mineral flux [g/m2*y]
r_grain          = 'PSD';      % mineral grain size | 'PSD' (Renforth distribution) | integer [um]
mineral          = 'basalt';   % mineral phase [string]
Calibration_site = 'GofM A';   % sites used for validation: Gulf of Mexico ('GofM A' , 'GofM B' ; https://doi.org/10.5194/bg-17-1685-2020), ...
%                                'Pacific'  (Growth and feeding of deep-sea coralLophelia pertusa from the
%                                Californiamargin under simulated ocean acidification conditions)
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- parameters for calibration ------------------------------------------
% -------------------------------------------------------------------------
if strcmp(Calibration_site,'GofM A')
% -------------------------------------------------------------------------    
F_NPP         = 1000;                    % net primary production NPP [gC/m2*y]
F_dust        = 2;                       % dust flux [g/m2*y]
mix_depth     = 50;                      % temperature mixing depth [m]
thermo_depth  = 80;                      % thermocline depth [m]
temperature_validation = temp_1_Validation;
DICinit       = 1000.*DIC_1_Validation(1,1);                    % dissolved inorganic carbon
TAinit        = 1000.*ALK_1_Validation(1,1);                    % alkalinity
Kz_surface    = 5E-3;                    % eddy diffusivity at surface [m2/s]
Kz_thermo     = 5E-5;                    % eddy diffusivity at thermocline [m2/s]
Kz_bottom     = 3E-4;                    % eddy diffusivity at bottom [m2/s]
K1_Adv_TH2    = 30;                      % fitting parameter for ALK advective flux
K2_Adv_TH2    = -4;                      % fitting parameter for ALK advective flux          
K_Adv_TH1     = 0.001;                   % fitting parameter for DIC advective flux
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------
if strcmp(Calibration_site,'GofM B')
% -------------------------------------------------------------------------
F_NPP         = 1000;                   % net primary production NPP [gC/m2*y]
F_dust        = 2;                      % dust flux [g/m2*y]
mix_depth     = 20;                     % temperature mixing depth [m]
thermo_depth  = 40;                     % thermocline depth [m]
temperature_validation = temp_2_Validation;
DICinit       = 1000.*DIC_2_Validation(1,1);                    % dissolved inorganic carbon
TAinit        = 1000.*ALK_2_Validation(1,1);                    % alkalinity
Kz_surface    = 5E-3;                   % eddy diffusivity at surface [m2/s]
Kz_thermo     = 5E-5;                   % eddy diffusivity at thermocline [m2/s]
Kz_bottom     = 6E-4;                   % eddy diffusivity at bottom [m2/s]
K1_Adv_TH2    = 38;                     % fitting parameter for ALK advective flux
K2_Adv_TH2    = -10;                    % fitting parameter for ALK advective flux          
K_Adv_TH1     = 0.00001;                % fitting parameter for DIC advective flux
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------
if strcmp(Calibration_site,'Pacific')
% -------------------------------------------------------------------------
F_NPP         = 1000;                   % net primary production NPP [gC/m2*y]
F_dust        = 1;                      % dust flux [g/m2*y]
mix_depth     = 20;                     % temperature mixing depth [m]
thermo_depth  = 40;                     % thermocline depth [m]
temperature_validation = Cali_temp;
DICinit       = 1800;                   % dissolved inorganic carbon
TAinit        = 0.999.*Cali_ALK(1,1);   % alkalinity
Kz_surface    = 5E-3;                   % eddy diffusivity at surface [m2/s]
Kz_thermo     = 4E-5;                   % eddy diffusivity at thermocline [m2/s]
Kz_bottom     = 5E-4;                   % eddy diffusivity at bottom [m2/s]
K1_Adv_TH2    = 0.025;                  % fitting parameter for ALK advective flux
K2_Adv_TH2    = 0.025;                  % fitting parameter for ALK advective flux          
K_Adv_TH1     = 0.00001;                % fitting parameter for DIC advective flux
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- fractional production terms -----------------------------------------
% -------------------------------------------------------------------------
K_pico        = 0.4;       % picoplankton fraction of primary production [dimensionless]
K_cocc        = 0.3;       % coccolith (eukaryotic) fraction of primary production [dimensionless]
K_arag        = 0.3;       % aragonite-containing phytoplankton fraction of primary production [dimensionless]
K_diat        = 0.0;       % diatom fraction of primary production [dimensionless]
K_dino        = 0.0;       % dinoflagellate fraction of primary production [dimensionless]
% ------------------------------------------------------------------------- 
if strcmp(r_grain,'PSD')
    % ---------------------------------------------------------------------
    % --- fraction of grains in each radius bin [dimensionless] -----------
    % ---------------------------------------------------------------------
    % NOTE: distribution here follows that of Renforth [2012] | doi:10.1016/j.ijggc.2012.06.011
    K_mineral_1  = 0.127; 
    K_mineral_2  = 0.139;
    K_mineral_3  = 0.108;
    K_mineral_4  = 0.102;    
    K_mineral_5  = 0.111;
    K_mineral_6  = 0.128; 
    K_mineral_7  = 0.062;    
    K_mineral_8  = 0.104;
    K_mineral_9  = 0.041;
    K_mineral_10 = 0.042;
    K_mineral_11 = 0.011;
    K_mineral_12 = 0.008;
    K_mineral_13 = 0.009;
    K_mineral_14 = 0.008;
    % ---------------------------------------------------------------------
    % --- grain radii [um] ------------------------------------------------
    % ---------------------------------------------------------------------
    r_mineral_1  =  0.09;
    r_mineral_2  =  0.10;
    r_mineral_3  =  0.50;
    r_mineral_4  =  1.00;
    r_mineral_5  =  2.00;
    r_mineral_6  =  4.00;
    r_mineral_7  =  8.63;
    r_mineral_8  = 12.20;
    r_mineral_9  = 17.30;    
    r_mineral_10 = 24.40;
    r_mineral_11 = 34.50;
    r_mineral_12 = 48.80;
    r_mineral_13 = 69.10;
    r_mineral_14 = 97.70;
else
    K_mineral = 1.0;
    r_mineral = r_grain;
end
% -------------------------------------------------------------------------
% --- misc parameters -----------------------------------------------------
% -------------------------------------------------------------------------
b             = 0.858;          % 'Martin' exponent for organic matter decay w/depth
% -------------------------------------------------------------------------
fractal_dim   = 2.0;            % aggregate fractal dimension
k_diss_arag   = 1E-3;          % rate constant for aragonite dissolution
k_diss_calc   = 1E-3;          % rate constant for calcite dissolution
num_small     = 1.e-20;         % arbitrarily small integer
% -------------------------------------------------------------------------
SST           = 25;             % sea surface temperature [degC]
T_bottom      = 5;              % temperature at ocean bottom [degC]
T_thermo      = 6;              % temperature at the bottom of thermocline
SSD           = 1.0255;         % sea surface density [g/cm3]
rho_bottom    = 1.0258;         % density at ocean bottom [g/cm3]
epsilon       = 1.e-4;          % turbulent dissipation rate [m2/s3]
mu            = 0.00188;        % kinematic viscosity of seawater [Ns/m2]
k_boltz       = 1.38064852e-23; % Boltzman constant [m2*kg/s2*K]
kernel_ref    = 5.e-18;         % reference kernel [m3/s]
F_NPP_ref     = 200;            % reference primary production [gC/m2/y]
F_dust_ref    = 10;             % reference dust flux [gC/m2/y]
Q10           = 2;              % temperature dependency factor
a_powerlaw    = 0.977;          % first power law coefficient (K&C 2015)
b_powerlaw    = 0.312;          % second power law coefficient (K&C 2015)
R             = 8.314;          % gas constant [J/K/mol]
C_to_K        = 273.15;         % temperature conversion - Celcius to Kelvin [1]
%
% -------------------------------------------------------------------------
% --- densities of primary particle components ----------------------------
% -------------------------------------------------------------------------
% NOTE: see Jokulsdottir + Archer [2016] | doi:10.5194/gmd-9-1455-2016
% NOTE: see Muscatine et al. [1998] | doi:10.1077/s003380050133
% NOTE: all in g/cm3
%
rho_org        = 1.06;          % fresh organic matter density
rho_carb_cocc  = 2.71;          % calcite density
rho_carb_arag  = 2.71;          % aragonite density
rho_bSi        = 2.09;          % biogenic silicate density
rho_dust       = 2.65;          % dust density
rho_dino       = 1.60;          % dinoflagellate density
%
% -------------------------------------------------------------------------
% --- silicate parameters -------------------------------------------------
% -------------------------------------------------------------------------
%
if strcmp(mineral,'basalt')
    rho_mineral      =  2.9;      % mineral density [g/cm3]
    k_mineral        =  1e-3;     % mineral dissolution rate [mol/m2*y]
    mineral_mass     =  215.0;    % molar mass of mineral [g/mol]
    mineral_porosity =  0.9;      % mineral porosity
    mineral_frac_dim =  2.0;      % mineral fractal dimension
    R_CO2            =  3.0;      % ideal capture stoichiometry [mol_CO2/mol_mineral]
    p1               = -0.02383;  % (-0.02536, -0.0223) coefficient for pH dependency (p1*x^3 + p2*x^2 + p3*x + p4; SSE: 0.64; R-square: 0.9874)
    p2               =  0.5625;   % (0.5333, 0.5916)
    p3               = -4.14;     % (-4.312, -3.969)
    p4               = -1.046;    % (-1.354, -0.7382)
    Temp_ref         =  25;       % reference temperature for dissolution rate constant
    A_0a             = 10^0.2;    % pre-exponential factor (acid) [mol/m2*s]
    A_0n             = 10^-4.3;   % pre-exponential factor (neutral) [mol/m2*s]
    A_0b             = 10^-0.45;  % pre-exponential factor (basic) [mol/m2*s]
    E_aa             = 32E3;      % apparent activation energy (acid) [KJ/mol]
    E_an             = 30E3;      % apparent activation energy (neutral) [KJ/mol]
    E_ab             = 30E3;      % apparent activation energy (basic) [KJ/mol]
    n_a              = 0.8;       % reaction order (acid) [1]
    n_n              = 0;         % reaction order (neutral) [1]
    n_b              = 0.7;       % reaction order (basic) [1]
end
% -------------------------------------------------------------------------
if strcmp(mineral,'olivine')
    rho_mineral      =  3.3;     % mineral density [g/cm3]
    k_mineral        =  3e-3;    % mineral dissolution rate [mol/m2*y]
    mineral_mass     =  214.0;   % molar mass of mineral [g/mol]
    mineral_porosity =  0.9;     % mineral porosity
    mineral_frac_dim =  2.0;     % mineral fractal dimension
    R_CO2            =  4.0;     % ideal capture stoichiometry [mol_CO2/mol_mineral]
    p1               = -0.3713;  %(-0.3744, -0.3682)
    p2               =  -7.271;  %(-7.293, -7.25)
    Temp_ref         =  25;      % reference temperature for dissolution rate constant
    A_0a             = 10^5.17;  % pre-exponential factor (acid) [mol/m2*s]
    A_0b             = 10^2.35;  % pre-exponential factor (basic) [mol/m2*s]
    E_aa             = 70.4E3;   % apparent activation energy (acid) [KJ/mol]
    E_ab             = 60.9E3;   % apparent activation energy (basic) [KJ/mol]
    n_a              = 0.46;     % reaction order (acid) [1]
    n_b              = 0.256;    % reaction order (basic) [1]
end
% -------------------------------------------------------------------------
if strcmp(mineral,'MgO')
    rho_mineral      = 3.6;      % mineral density [g/cm3]
    k_mineral        = 4e-1;     % mineral dissolution rate [mol/m2*y]
    mineral_mass     = 40.0;     % molar mass of mineral [g/mol]
    mineral_porosity = 0.9;      % mineral porosity
    mineral_frac_dim = 2.0;      % mineral fractal dimension
    R_CO2            = 1.0;      % ideal capture stoichiometry [mol_CO2/mol_mineral]
    Temp_ref         = 25;       % reference temperature for dissolution rate constant
    A_0a             = 10^3.57;  % pre-exponential factor (acid) [mol/m2*s]
    A_0n             = -10^-0.77;% pre-exponential factor (neutral) [mol/m2*s]
    E_aa             = 59E3;     % apparent activation energy (acid) [KJ/mol]
    E_an             = 42E3;     % apparent activation energy (neutral) [KJ/mol]
    n_a              = 0.115;    % reaction order (acid) [1]
    n_n              = 0;        % reaction order (neutral) [1]
end
% -------------------------------------------------------------------------
if strcmp(mineral,'CaO')
    rho_mineral      = 3.3;      % mineral density [g/cm3]
    k_mineral        = 1.6;      % mineral dissolution rate [mol/m2*y]
    mineral_mass     = 56.1;     % molar mass of mineral [g/mol]
    mineral_porosity = 0.9;      % mineral porosity
    mineral_frac_dim = 2.0;      % mineral fractal dimension
    R_CO2            = 1.0;      % ideal capture stoichiometry [mol_CO2/mol_mineral]
    Temp_ref         = 25;       % reference temperature for dissolution rate constant
    A_0              = 10^-0.72; % pre-exponential factor [mol/m2*s]
    E_a              = 13.72E3;  % apparent activation energy [KJ/mol]
    n_all            = 0.0438;   % reaction order [1]
end
%
% -------------------------------------------------------------------------
% --- calculate SSA and dissolution rate constant values ------------------
% -------------------------------------------------------------------------
% NOTE: SSA in units of m2/g
% NOTE: SSA calculated based on grain size according to original arrays
if strcmp(r_grain,'PSD')
    SSA_1             = (15.53*(r_mineral_1^-0.2829))  - 1.437;
    SSA_2             = (15.53*(r_mineral_2^-0.2829))  - 1.437;
    SSA_3             = (15.53*(r_mineral_3^-0.2829))  - 1.437;
    SSA_4             = (15.53*(r_mineral_4^-0.2829))  - 1.437;
    SSA_5             = (15.53*(r_mineral_5^-0.2829))  - 1.437;
    SSA_6             = (15.53*(r_mineral_6^-0.2829))  - 1.437;
    SSA_7             = (15.53*(r_mineral_7^-0.2829))  - 1.437;
    SSA_8             = (15.53*(r_mineral_8^-0.2829))  - 1.437;
    SSA_9             = (15.53*(r_mineral_9^-0.2829))  - 1.437;
    SSA_10            = (15.53*(r_mineral_10^-0.2829)) - 1.437;
    SSA_11            = (15.53*(r_mineral_11^-0.2829)) - 1.437;
    SSA_12            = (15.53*(r_mineral_12^-0.2829)) - 1.437;
    SSA_13            = (15.53*(r_mineral_13^-0.2829)) - 1.437;
    SSA_14            = (15.53*(r_mineral_14^-0.2829)) - 1.437;
else
    SSA               = (15.53*(r_mineral^-0.2829)) - 1.437;
end
% -------------------------------------------------------------------------
% NOTE: rate constants in units of /y
if strcmp(r_grain,'PSD')
    k_diss_mineral_1  = k_mineral.*SSA_1.*mineral_mass;
    k_diss_mineral_2  = k_mineral.*SSA_2.*mineral_mass;
    k_diss_mineral_3  = k_mineral.*SSA_3.*mineral_mass;
    k_diss_mineral_4  = k_mineral.*SSA_4.*mineral_mass;
    k_diss_mineral_5  = k_mineral.*SSA_5.*mineral_mass;
    k_diss_mineral_6  = k_mineral.*SSA_6.*mineral_mass;
    k_diss_mineral_7  = k_mineral.*SSA_7.*mineral_mass;
    k_diss_mineral_8  = k_mineral.*SSA_8.*mineral_mass;
    k_diss_mineral_9  = k_mineral.*SSA_9.*mineral_mass;
    k_diss_mineral_10 = k_mineral.*SSA_10.*mineral_mass;
    k_diss_mineral_11 = k_mineral.*SSA_11.*mineral_mass;
    k_diss_mineral_12 = k_mineral.*SSA_12.*mineral_mass;
    k_diss_mineral_13 = k_mineral.*SSA_13.*mineral_mass;
    k_diss_mineral_14 = k_mineral.*SSA_14.*mineral_mass;
else
    k_diss_mineral    = k_mineral.*SSA.*mineral_mass;
end
%
% -------------------------------------------------------------------------
% --- calculate characteristics of particle size classes ------------------
% -------------------------------------------------------------------------
% NOTE: correct particle radii for fractal dimension
% NOTE: radius values in cm
if strcmp(r_grain,'PSD')
    r_mineral_1   = (r_mineral_1*1.e-4)^(3/fractal_dim);
    r_mineral_2   = (r_mineral_2*1.e-4)^(3/fractal_dim);
    r_mineral_3   = (r_mineral_3*1.e-4)^(3/fractal_dim);
    r_mineral_4   = (r_mineral_4*1.e-4)^(3/fractal_dim);
    r_mineral_5   = (r_mineral_5*1.e-4)^(3/fractal_dim);
    r_mineral_6   = (r_mineral_6*1.e-4)^(3/fractal_dim);
    r_mineral_7   = (r_mineral_7*1.e-4)^(3/fractal_dim);
    r_mineral_8   = (r_mineral_8*1.e-4)^(3/fractal_dim);
    r_mineral_9   = (r_mineral_9*1.e-4)^(3/fractal_dim);
    r_mineral_10  = (r_mineral_10*1.e-4)^(3/fractal_dim);
    r_mineral_11  = (r_mineral_11*1.e-4)^(3/fractal_dim);
    r_mineral_12  = (r_mineral_12*1.e-4)^(3/fractal_dim);
    r_mineral_13  = (r_mineral_13*1.e-4)^(3/fractal_dim);
    r_mineral_14  = (r_mineral_14*1.e-4)^(3/fractal_dim);
else
    r_mineral     = (r_mineral*1.e-4)^(3/fractal_dim);
end
% -------------------------------------------------------------------------
% NOTE: all mineral volumes in cm3
if strcmp(r_grain,'PSD')
    v_mineral_1   = (4/3)*pi*(r_mineral_1^3);
    v_mineral_2   = (4/3)*pi*(r_mineral_2^3);  
    v_mineral_3   = (4/3)*pi*(r_mineral_3^3);  
    v_mineral_4   = (4/3)*pi*(r_mineral_4^3);  
    v_mineral_5   = (4/3)*pi*(r_mineral_5^3);  
    v_mineral_6   = (4/3)*pi*(r_mineral_6^3);  
    v_mineral_7   = (4/3)*pi*(r_mineral_7^3); 
    v_mineral_8   = (4/3)*pi*(r_mineral_8^3); 
    v_mineral_9   = (4/3)*pi*(r_mineral_9^3); 
    v_mineral_10  = (4/3)*pi*(r_mineral_10^3); 
    v_mineral_11  = (4/3)*pi*(r_mineral_11^3); 
    v_mineral_12  = (4/3)*pi*(r_mineral_12^3); 
    v_mineral_13  = (4/3)*pi*(r_mineral_13^3); 
    v_mineral_14  = (4/3)*pi*(r_mineral_14^3); 
else
    v_mineral     = (4/3)*pi*(r_mineral^3);
end
% -------------------------------------------------------------------------
% NOTE: all mineral masses in g/particle
if strcmp(r_grain,'PSD')
    m_mineral_1   = v_mineral_1  * rho_mineral;
    m_mineral_2   = v_mineral_2  * rho_mineral;
    m_mineral_3   = v_mineral_3  * rho_mineral;
    m_mineral_4   = v_mineral_4  * rho_mineral;
    m_mineral_5   = v_mineral_5  * rho_mineral;
    m_mineral_6   = v_mineral_6  * rho_mineral;
    m_mineral_7   = v_mineral_7  * rho_mineral;
    m_mineral_8   = v_mineral_8  * rho_mineral;
    m_mineral_9   = v_mineral_9  * rho_mineral;
    m_mineral_10  = v_mineral_10 * rho_mineral;
    m_mineral_11  = v_mineral_11 * rho_mineral;
    m_mineral_12  = v_mineral_12 * rho_mineral;
    m_mineral_13  = v_mineral_13 * rho_mineral;
    m_mineral_14  = v_mineral_14 * rho_mineral;
else
    m_mineral     = v_mineral    * rho_mineral;
end
%
% -------------------------------------------------------------------------
% --- primary particle masses ---------------------------------------------
% -------------------------------------------------------------------------
% NOTE: see Jokulsdottir + Archer [2016] | doi:10.5194/gmd-9-1455-2016
% NOTE: all in g/particle
% -------------------------------------------------------------------------
m_org_pico    = 1.00e-12;       % organic mass in picoplankton particle
% -------------------------------------------------------------------------
m_org_cocc    = 7.00e-12;       % organic mass in coccolithophore particle
m_carb_cocc   = 1.50e-12;       % calcite mass in coccolithophore particle
% -------------------------------------------------------------------------
m_org_arag    = 7.00e-12;       % organic mass in aragonite-producing phytoplankton particle
m_carb_arag   = 1.50e-12;       % aragonite mass in aragonite-producing phytoplankton particle
% -------------------------------------------------------------------------
m_org_diat    = 15.0e-12;       % organic mass in diatom particle
m_bSi_diat    = 5.00e-12;       % biogenic silica mass in diatom particle
% -------------------------------------------------------------------------
m_org_dino    = 4.00e-12;       % organic mass in dinoflagellate particle 
% -------------------------------------------------------------------------
m_dust        = 2.85e-12;       % mass of dust in primary particle
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- INITIALIZATION ------------------------------------------------------
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- set up grids for depth | density | temperature | viscosity ----------
% -------------------------------------------------------------------------
n_depth       = 500;                                % number of depth intervals
max_depth     = 500;                                % max ocean depth [m]
dz            = max_depth/n_depth;                  % ocean depth interval [m]
z             = linspace(1,max_depth,n_depth);      % ocean depth grid [m]
% -------------------------------------------------------------------------
rho_seawater                                     = ones(1,n_depth);
rho_seawater(1,1:mix_depth/(max_depth/n_depth))                   = SSD;
rho_seawater(1,(mix_depth/(max_depth/n_depth)+1):thermo_depth/(max_depth/n_depth)) = ( (rho_bottom - SSD)./ ...
                                                    (thermo_depth-mix_depth) ).* ...
                                                   z(1,(mix_depth/(max_depth/n_depth)+1):thermo_depth/(max_depth/n_depth))+ ...
                                                   SSD- ...
                                                   ((rho_bottom - SSD)./(thermo_depth-mix_depth)).* ...
                                                   z(1,(mix_depth/(max_depth/n_depth)+1));
rho_seawater(1,(thermo_depth/(max_depth/n_depth)+1):n_depth)          = rho_bottom*ones(1,n_depth-thermo_depth/(max_depth/n_depth));
% -------------------------------------------------------------------------
T_seawater = interp1(temperature_validation(:,2),temperature_validation(:,1),z); 
% -------------------------------------------------------------------------
Depth_scale = thermo_depth;                         % thermocline depth      
K_temperature             = (Q10.^((T_seawater-Temp_ref)./(max_depth/n_depth)));   % depth profile of temperature dependency factor
% -------------------------------------------------------------------------
nho                       = 5.e-5 * exp(2250./(T_seawater+273.15));
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- initial primary particle densities for 'mixed' particles ------------
% -------------------------------------------------------------------------
f_org_cocc    = m_org_cocc/(m_org_cocc+m_carb_cocc);
f_carb_cocc   = m_carb_cocc/(m_org_cocc+m_carb_cocc);
% -------------------------------------------------------------------------
f_org_arag    = m_org_arag/(m_org_arag+m_carb_arag);
f_carb_arag   = m_carb_arag/(m_org_arag+m_carb_arag);
% -------------------------------------------------------------------------
f_org_diat    = m_org_diat/(m_org_diat+m_bSi_diat);
f_bSi_diat    = m_bSi_diat/(m_org_diat+m_bSi_diat);
%
rho_pico      = rho_org;
rho_cocc      = (f_org_cocc*rho_org)+(f_carb_cocc*rho_carb_cocc);
rho_arag      = (f_org_arag*rho_org)+(f_carb_arag*rho_carb_arag);
rho_diat      = (f_org_diat*rho_org)+(f_bSi_diat*rho_bSi);
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- calculate material fluxes -------------------------------------------
% -------------------------------------------------------------------------
Flux_pico        = K_pico.*F_NPP;                   % picoplankton flux [gC/m2*y]
Flux_org_cocc    = f_org_cocc.*K_cocc.*F_NPP;       % coccolith organic flux [gC/m2*y]
Flux_carb_cocc   = f_carb_cocc.*K_cocc.*F_NPP;      % coccolith carbonate (calcite) flux [g/m2*y]
Flux_org_arag    = f_org_arag.*K_arag.*F_NPP;       % aragonite organic flux [g/m2*y]
Flux_carb_arag   = f_carb_arag.*K_arag.*F_NPP;      % aragonite carbonate (aragonite) flux [g/m2*y]
Flux_org_diat    = f_org_diat.*K_diat.*F_NPP;       % diatom organic flux [gC/m2*y]
Flux_bSi_diat    = f_bSi_diat.*K_diat.*F_NPP;       % diatom biogenic Si flux
Flux_dino        = K_dino.*F_NPP;                   % dino flux [gC/m2*y]
% -------------------------------------------------------------------------
if strcmp(r_grain,'PSD')
    Flux_mineral_1   = K_mineral_1.*F_mineral;
    Flux_mineral_2   = K_mineral_2.*F_mineral;
    Flux_mineral_3   = K_mineral_3.*F_mineral;
    Flux_mineral_4   = K_mineral_4.*F_mineral;
    Flux_mineral_5   = K_mineral_5.*F_mineral;
    Flux_mineral_6   = K_mineral_6.*F_mineral;
    Flux_mineral_7   = K_mineral_7.*F_mineral;
    Flux_mineral_8   = K_mineral_8.*F_mineral;
    Flux_mineral_9   = K_mineral_9.*F_mineral;
    Flux_mineral_10  = K_mineral_10.*F_mineral;
    Flux_mineral_11  = K_mineral_11.*F_mineral;
    Flux_mineral_12  = K_mineral_12.*F_mineral;
    Flux_mineral_13  = K_mineral_13.*F_mineral;
    Flux_mineral_14  = K_mineral_14.*F_mineral;
    Flux_mineral_tot =  Flux_mineral_1  + Flux_mineral_2  + Flux_mineral_3 ...
                      + Flux_mineral_4  + Flux_mineral_5  + Flux_mineral_6 ...
                      + Flux_mineral_7  + Flux_mineral_8  + Flux_mineral_9 ...
                      + Flux_mineral_10 + Flux_mineral_11 + Flux_mineral_12 ...
                      + Flux_mineral_13 + Flux_mineral_14;
else
    Flux_mineral     = K_mineral*F_mineral;
    Flux_mineral_tot = Flux_mineral;
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% --- initisl depth profiles of material fluxes ---------------------------
% -------------------------------------------------------------------------
Flux_pico_z = Flux_pico.*(z.^-b);
Flux_org_cocc_z = Flux_org_cocc.*(z.^-b);
Flux_carb_cocc_z = max(Flux_carb_cocc.*-k_diss_calc.*z.*1.e-10,num_small);
Flux_org_arag_z = Flux_org_arag.*(z.^-b);
Flux_carb_arag_z = max(Flux_carb_arag.*-k_diss_calc.*z.*1.e-10,num_small);
Flux_org_diat_z = Flux_org_diat.*(z.^-b);
Flux_bSi_diat_z = ones(1,n_depth).*Flux_bSi_diat;
Flux_dust_z = ones(1,n_depth).*F_dust;
Flux_dino_z = Flux_dino.*(z.^-b);
% -------------------------------------------------------------------------
if strcmp(r_grain,'PSD')
    Flux_mineral_z_1  = max(Flux_mineral_1  - (k_diss_mineral_1.*z*1.e-10),num_small); 
    Flux_mineral_z_2  = max(Flux_mineral_2  - (k_diss_mineral_2.*z*1.e-10),num_small);
    Flux_mineral_z_3  = max(Flux_mineral_3  - (k_diss_mineral_3.*z*1.e-10),num_small);
    Flux_mineral_z_4  = max(Flux_mineral_4  - (k_diss_mineral_4.*z*1.e-10),num_small);
    Flux_mineral_z_5  = max(Flux_mineral_5  - (k_diss_mineral_5.*z*1.e-10),num_small);
    Flux_mineral_z_6  = max(Flux_mineral_6  - (k_diss_mineral_6.*z*1.e-10),num_small);
    Flux_mineral_z_7  = max(Flux_mineral_7  - (k_diss_mineral_7.*z*1.e-10),num_small);
    Flux_mineral_z_8  = max(Flux_mineral_8  - (k_diss_mineral_8.*z*1.e-10),num_small);
    Flux_mineral_z_9  = max(Flux_mineral_9  - (k_diss_mineral_9.*z*1.e-10),num_small);
    Flux_mineral_z_10 = max(Flux_mineral_10 - (k_diss_mineral_10.*z*1.e-10),num_small);
    Flux_mineral_z_11 = max(Flux_mineral_11 - (k_diss_mineral_11.*z*1.e-10),num_small);
    Flux_mineral_z_12 = max(Flux_mineral_12 - (k_diss_mineral_12.*z*1.e-10),num_small);
    Flux_mineral_z_13 = max(Flux_mineral_13 - (k_diss_mineral_13.*z*1.e-10),num_small);
    Flux_mineral_z_14 = max(Flux_mineral_14 - (k_diss_mineral_14.*z*1.e-10),num_small);
    Flux_mineral_tot_z =  Flux_mineral_z_1  + Flux_mineral_z_2  + Flux_mineral_z_3 ...
                          + Flux_mineral_z_4  + Flux_mineral_z_5  + Flux_mineral_z_6 ...
                          + Flux_mineral_z_7  + Flux_mineral_z_8  + Flux_mineral_z_9 ...
                          + Flux_mineral_z_10 + Flux_mineral_z_11 + Flux_mineral_z_12 ...
                          + Flux_mineral_z_13 + Flux_mineral_z_14;
else
    Flux_mineral_z    = max(Flux_mineral    - (k_diss_mineral*z*1.e-12),num_small);
    Flux_mineral_tot_z = Flux_mineral_z;
end

% -------------------------------------------------------------------------

Flux_total_z       =  Flux_pico_z ...
                      + Flux_org_cocc_z ...
                      + Flux_carb_cocc_z ...
                      + Flux_org_arag_z ...
                      + Flux_carb_arag_z ...
                      + Flux_org_diat_z ...
                      + Flux_bSi_diat_z ...
                      + Flux_dino_z ...
                      + Flux_dust_z ...
                      + Flux_mineral_tot_z;                            % total flux
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- calculate flux ratios for aggregation probability -------------------
% -------------------------------------------------------------------------
R_pico_z        = Flux_pico_z./Flux_total_z;
R_dino_z        = Flux_dino_z./Flux_total_z;
% -------------------------------------------------------------------------
R_org_cocc_z    = Flux_org_cocc_z./Flux_total_z;
R_carb_cocc_z   = Flux_carb_cocc_z./Flux_total_z;
R_cocc_z        = R_org_cocc_z+R_carb_cocc_z;
% -------------------------------------------------------------------------
R_org_arag_z    = Flux_org_arag_z./Flux_total_z;
R_carb_arag_z   = Flux_carb_arag_z./Flux_total_z;
R_arag_z        = R_org_arag_z+R_carb_arag_z;
% -------------------------------------------------------------------------
R_org_diat_z    = Flux_org_diat_z./Flux_total_z;
R_bSi_diat_z    = Flux_bSi_diat_z./Flux_total_z;
R_diat_z        = R_org_diat_z+R_bSi_diat_z;
% -------------------------------------------------------------------------
R_dust_z        = Flux_dust_z./Flux_total_z;
% -------------------------------------------------------------------------
if strcmp(r_grain,'PSD')
    R_mineral_z_1   = Flux_mineral_z_1./Flux_total_z;
    R_mineral_z_2   = Flux_mineral_z_2./Flux_total_z;
    R_mineral_z_3   = Flux_mineral_z_3./Flux_total_z;
    R_mineral_z_4   = Flux_mineral_z_4./Flux_total_z;
    R_mineral_z_5   = Flux_mineral_z_5./Flux_total_z;
    R_mineral_z_6   = Flux_mineral_z_6./Flux_total_z;
    R_mineral_z_7   = Flux_mineral_z_7./Flux_total_z;
    R_mineral_z_8   = Flux_mineral_z_8./Flux_total_z;
    R_mineral_z_9   = Flux_mineral_z_9./Flux_total_z;
    R_mineral_z_10  = Flux_mineral_z_10./Flux_total_z;
    R_mineral_z_11  = Flux_mineral_z_11./Flux_total_z;
    R_mineral_z_12  = Flux_mineral_z_12./Flux_total_z;
    R_mineral_z_13  = Flux_mineral_z_13./Flux_total_z;
    R_mineral_z_14  = Flux_mineral_z_14./Flux_total_z;
else
    R_mineral_z     = Flux_mineral_z./Flux_total_z;
end
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- initial mass calculations for materials -----------------------------
% -------------------------------------------------------------------------
m_org_pico_z  = m_org_pico*(z.^-b);                            % picoplankton organic mass at depth [g/particle]
m_org_cocc_z  = m_org_cocc*(z.^-b);                            % coccolith organic mass at depth [g/particle]
m_org_arag_z  = m_org_arag*(z.^-b);                            % aragonite-producing phytoplankton organic mass at depth [g/particle]
m_org_diat_z  = m_org_diat*(z.^-b);                            % diatom organic mass at depth [g/particle]
m_org_dino_z  = m_org_dino*(z.^-b);                            % dinoflagellate organic mass at depth [g/particle]
% -------------------------------------------------------------------------
m_carb_cocc_z = max(m_carb_cocc-k_diss_calc.*z.*1.e-12,num_small); % coccolith calcite mass at depth [g/particle]
m_carb_arag_z = max(m_carb_arag-k_diss_arag.*z*1.e-12,num_small);  % aragonite-producing phytoplankton aragonite mass at depth [g/particle]
m_bSi_diat_z  = ones(1,n_depth).*m_bSi_diat;                   % diatom biogenic silica mass at depth [g/particle]
% -------------------------------------------------------------------------
F_NPP_z = F_NPP*(z.^-b);
F_NPP_ref_z = F_NPP_ref*(z.^-b);
% -------------------------------------------------------------------------
% NOTE: all silicate mineral masses in g/particle
if strcmp(r_grain,'PSD')
    m_mineral_1_z  = max(m_mineral_1  - (k_diss_mineral_1.*z*1.e-12),num_small); 
    m_mineral_2_z  = max(m_mineral_2  - (k_diss_mineral_2.*z*1.e-12),num_small);
    m_mineral_3_z  = max(m_mineral_3  - (k_diss_mineral_3.*z*1.e-12),num_small);
    m_mineral_4_z  = max(m_mineral_4  - (k_diss_mineral_4.*z*1.e-12),num_small);
    m_mineral_5_z  = max(m_mineral_5  - (k_diss_mineral_5.*z*1.e-12),num_small);
    m_mineral_6_z  = max(m_mineral_6  - (k_diss_mineral_6.*z*1.e-12),num_small);
    m_mineral_7_z  = max(m_mineral_7  - (k_diss_mineral_7.*z*1.e-12),num_small);
    m_mineral_8_z  = max(m_mineral_8  - (k_diss_mineral_8.*z*1.e-12),num_small);
    m_mineral_9_z  = max(m_mineral_9  - (k_diss_mineral_9.*z*1.e-12),num_small);
    m_mineral_10_z = max(m_mineral_10 - (k_diss_mineral_10.*z*1.e-12),num_small);
    m_mineral_11_z = max(m_mineral_11 - (k_diss_mineral_11.*z*1.e-12),num_small);
    m_mineral_12_z = max(m_mineral_12 - (k_diss_mineral_12.*z*1.e-12),num_small);
    m_mineral_13_z = max(m_mineral_13 - (k_diss_mineral_13.*z*1.e-12),num_small);
    m_mineral_14_z = max(m_mineral_14 - (k_diss_mineral_14.*z*1.e-12),num_small);
else
    m_mineral_z    = max(m_mineral    - (k_diss_mineral*z*1.e-12),num_small);
end
%
% -------------------------------------------------------------------------              
%              
% -------------------------------------------------------------------------
% --- determining grain size of most probable size class
% -------------------------------------------------------------------------
if strcmp(r_grain,'PSD')
    max_mineral = max([K_mineral_1  K_mineral_2  K_mineral_3  K_mineral_4 ...
                       K_mineral_5  K_mineral_6  K_mineral_7  K_mineral_8 ...
                       K_mineral_9  K_mineral_10 K_mineral_11 K_mineral_12 ...
                       K_mineral_13 K_mineral_14]);
    % ---------------------------------------------------------------------
    if K_mineral_1  == max_mineral
        m_mineral_max_z = m_mineral_1_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_2  == max_mineral
        m_mineral_max_z = m_mineral_2_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_3  == max_mineral
        m_mineral_max_z = m_mineral_3_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_4  == max_mineral
        m_mineral_max_z = m_mineral_4_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_5  == max_mineral
        m_mineral_max_z = m_mineral_5_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_6  == max_mineral
        m_mineral_max_z = m_mineral_6_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_7  == max_mineral
        m_mineral_max_z = m_mineral_7_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_8  == max_mineral
        m_mineral_max_z = m_mineral_8_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_9  == max_mineral
        m_mineral_max_z = m_mineral_9_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_10 == max_mineral
        m_mineral_max_z = m_mineral_10_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_11 == max_mineral
        m_mineral_max_z = m_mineral_11_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_12 ==  max_mineral
        m_mineral_max_z = m_mineral_12_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_13 == max_mineral
        m_mineral_max_z = m_mineral_13_z;
    end
    % ---------------------------------------------------------------------
    if K_mineral_14 == max_mineral
        m_mineral_max_z = m_mineral_14_z;
    end
else
    max_mineral         = K_mineral;
    m_mineral_max_z     = m_mineral_z;
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
m_pico_z            = m_org_pico_z;                  % total picoplankton mass
m_dino_z            = m_org_dino_z;                  % total dinoflagellate mass
m_cocc_tot_z        = m_org_cocc_z + m_carb_cocc_z;  % total coccolith mass
m_arag_tot_z        = m_org_arag_z + m_carb_arag_z;  % total aragonite plankton mass
m_diat_tot_z        = m_org_diat_z + m_bSi_diat_z;   % total diatom mass
% -------------------------------------------------------------------------
m_dust_tot_z        = ones(1,n_depth).*m_dust;           % total dust mass
% -------------------------------------------------------------------------
% NOTE: mass of silicate in each size class
if strcmp(r_grain,'PSD')
    m_mineral_tot_1_z   = m_mineral_1_z;
    m_mineral_tot_2_z   = m_mineral_2_z;
    m_mineral_tot_3_z   = m_mineral_3_z;
    m_mineral_tot_4_z   = m_mineral_4_z;
    m_mineral_tot_5_z   = m_mineral_5_z;
    m_mineral_tot_6_z   = m_mineral_6_z;
    m_mineral_tot_7_z   = m_mineral_7_z;
    m_mineral_tot_8_z   = m_mineral_8_z;
    m_mineral_tot_9_z   = m_mineral_9_z;
    m_mineral_tot_10_z  = m_mineral_10_z;
    m_mineral_tot_11_z  = m_mineral_11_z;
    m_mineral_tot_12_z  = m_mineral_12_z;
    m_mineral_tot_13_z  = m_mineral_13_z;
    m_mineral_tot_14_z  = m_mineral_14_z;
    m_mineral_tot_max_z = m_mineral_max_z;
else
    m_mineral_tot_z     = m_mineral_z;
    m_mineral_tot_max_z = m_mineral_max_z;
end
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- initialize volumes and densities for particle materials -------------
% -------------------------------------------------------------------------
% NOTE: fractions of each material [org|carb|bSi] in plankton
f_org_cocc_z     = m_org_cocc_z/(m_org_cocc_z+m_carb_cocc_z);
f_carb_cocc_z    = m_carb_cocc_z/(m_org_cocc_z+m_carb_cocc_z);
f_org_arag_z     = m_org_arag_z/(m_org_arag_z+m_carb_arag_z);
f_carb_arag_z    = m_carb_arag_z/(m_org_arag_z+m_carb_arag_z);
f_org_diat_z     = m_org_diat_z/(m_org_diat_z+m_bSi_diat_z);
f_bSi_diat_z     = m_bSi_diat_z/(m_org_diat_z+m_bSi_diat_z);
% -------------------------------------------------------------------------
% NOTE: densities of plankton calculated via mass balance
rho_pico_z       = rho_org;
rho_dino_z       = rho_dino;
rho_cocc_z       = (f_org_cocc_z*rho_org)+(f_carb_cocc_z*rho_carb_cocc);
rho_arag_z       = (f_org_arag_z*rho_org)+(f_carb_arag_z*rho_carb_arag);
rho_diat_z       = (f_org_diat_z*rho_org)+(f_bSi_diat_z*rho_bSi);
% -------------------------------------------------------------------------
rho_dust_z       = rho_dust;
% -------------------------------------------------------------------------
% NOTE: densities of all size classes defined by mineral phase
if strcmp(r_grain,'PSD')
    rho_mineral_1_z  = rho_mineral;
    rho_mineral_2_z  = rho_mineral;
    rho_mineral_3_z  = rho_mineral;
    rho_mineral_4_z  = rho_mineral;
    rho_mineral_5_z  = rho_mineral;
    rho_mineral_6_z  = rho_mineral;
    rho_mineral_7_z  = rho_mineral;
    rho_mineral_8_z  = rho_mineral;
    rho_mineral_9_z  = rho_mineral;
    rho_mineral_10_z = rho_mineral;
    rho_mineral_11_z = rho_mineral;
    rho_mineral_12_z = rho_mineral;
    rho_mineral_13_z = rho_mineral;
    rho_mineral_14_z = rho_mineral;
else
    rho_mineral_z    = rho_mineral;
end
% -------------------------------------------------------------------------
% NOTE: volumes of particle constituents [cm3]
v_pico_z        = m_pico_z/rho_pico_z;
v_cocc_z        = m_cocc_tot_z/rho_cocc_z;
v_arag_z        = m_arag_tot_z/rho_arag_z;
v_diat_z        = m_diat_tot_z/rho_diat_z;
v_dino_z        = m_dino_z/rho_dino_z;
% -------------------------------------------------------------------------
v_dust_z        = m_dust_tot_z/rho_dust_z;
% -------------------------------------------------------------------------
if strcmp(r_grain,'PSD')
    v_mineral_1_z   = m_mineral_tot_1_z/rho_mineral_1_z;
    v_mineral_2_z   = m_mineral_tot_2_z/rho_mineral_2_z;
    v_mineral_3_z   = m_mineral_tot_3_z/rho_mineral_3_z;
    v_mineral_4_z   = m_mineral_tot_4_z/rho_mineral_4_z;
    v_mineral_5_z   = m_mineral_tot_5_z/rho_mineral_5_z;
    v_mineral_6_z   = m_mineral_tot_6_z/rho_mineral_6_z;
    v_mineral_7_z   = m_mineral_tot_7_z/rho_mineral_7_z;
    v_mineral_8_z   = m_mineral_tot_8_z/rho_mineral_8_z;
    v_mineral_9_z   = m_mineral_tot_9_z/rho_mineral_9_z;
    v_mineral_10_z  = m_mineral_tot_10_z/rho_mineral_10_z;
    v_mineral_11_z  = m_mineral_tot_11_z/rho_mineral_11_z;
    v_mineral_12_z  = m_mineral_tot_12_z/rho_mineral_12_z;
    v_mineral_13_z  = m_mineral_tot_13_z/rho_mineral_13_z;
    v_mineral_14_z  = m_mineral_tot_14_z/rho_mineral_14_z;
    v_mineral_max_z = m_mineral_tot_max_z/rho_mineral_14_z;
else
    v_mineral_z     = m_mineral_tot_z/rho_mineral_z;
    v_mineral_max_z = m_mineral_tot_max_z/rho_mineral_z;
end
 
% -------------------------------------------------------------------------
for i=1:500
    
    i_1(1,i)             = round(R_pico_z(1,i).*4.e2);
    i_2(1,i)             = round(R_cocc_z(1,i).*4.e2);
    i_3(1,i)             = round(R_arag_z(1,i).*4.e2);
    i_4(1,i)             = round(R_diat_z(1,i).*4.e2);
    i_5(1,i)             = round(R_dust_z(1,i).*4.e2);
    i_6(1,i)             = round(R_dino_z(1,i).*4.e2);
    % -------------------------------------------------------------------------
    if strcmp(r_grain,'PSD')
        i_7(1,i)             = round(R_mineral_z_1(1,i).*4.e2);
        i_8(1,i)             = round(R_mineral_z_2(1,i).*4.e2);
        i_9(1,i)             = round(R_mineral_z_3(1,i).*4.e2);
        i_10(1,i)            = round(R_mineral_z_4(1,i).*4.e2);
        i_11(1,i)            = round(R_mineral_z_5(1,i).*4.e2);
        i_12(1,i)            = round(R_mineral_z_6(1,i).*4.e2);
        i_13(1,i)            = round(R_mineral_z_7(1,i).*4.e2);
        i_14(1,i)            = round(R_mineral_z_8(1,i).*4.e2);
        i_15(1,i)            = round(R_mineral_z_9(1,i).*4.e2);
        i_16(1,i)            = round(R_mineral_z_10(1,i).*4.e2);
        i_17(1,i)            = round(R_mineral_z_11(1,i).*4.e2);
        i_18(1,i)            = round(R_mineral_z_12(1,i).*4.e2);
        i_19(1,i)            = round(R_mineral_z_13(1,i).*4.e2);
        i_20(1,i)            = round(R_mineral_z_14(1,i).*4.e2);
    else
        i_7(1,i)             = round(R_mineral_z(1,i).*4.e2);
    end

end
%
for i_main=1:500
% -------------------------------------------------------------------------
    for i = 1:i_1(1,i_main)
        rho_particle(i_main,i) = rho_pico_z;
    end
    % ---------------------------------------------------------------------
    for j = ( i_1(1,i_main) + 1 ):...
            ( i_2(1,i_main) + i_1(1,i_main) + 1 )
        rho_particle(i_main,j) = rho_cocc_z;
    end
    % ---------------------------------------------------------------------
    for k = ( i_2(1,i_main) + i_1(1,i_main) + 2 ):...
            ( i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 2 )
        rho_particle(i_main,k) = rho_arag_z;
    end
    % ---------------------------------------------------------------------
    for m = ( i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 3 ):...
            ( i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 3 )
        rho_particle(i_main,m) = rho_diat_z;
    end
    % ---------------------------------------------------------------------
    for q = ( i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 4 ):...
            ( i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 4 )
        rho_particle(i_main,q) = rho_dust_z;
    end
    % ---------------------------------------------------------------------
    for pp = ( i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 5 ):...
             ( i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 5 )
        rho_particle(i_main,pp) = rho_dino_z;
    end
    % ---------------------------------------------------------------------
    if strcmp(r_grain,'PSD')
        for Min1 = ( i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                    +i_1(1,i_main) + 6 ):...
                   ( i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                    +i_2(1,i_main) + i_1(1,i_main) + 6 )
            rho_particle(i_main,Min1) = rho_mineral_1_z;
        end
        % -----------------------------------------------------------------
        for Min2 = ( i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                    +i_2(1,i_main) + i_1(1,i_main) + 7 ):...
                   ( i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                    +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 7 )
            rho_particle(i_main,Min2) = rho_mineral_2_z;
        end
        % -----------------------------------------------------------------
        for Min3 = ( i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                    +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 8 ):...
                   ( i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                    +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 8 )
            rho_particle(i_main,Min3) = rho_mineral_3_z;
        end
        % -----------------------------------------------------------------
        for Min4 = ( i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                    +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 9 ):...
                   ( i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main)...
                    +i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 9)
            rho_particle(i_main,Min4) = rho_mineral_4_z;
        end
        % -----------------------------------------------------------------
        for Min5 = ( i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main)...
                    +i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 10):...
                   ( i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main)...
                    +i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                    +i_1(1,i_main) + 10 )
            rho_particle(i_main,Min5) = rho_mineral_5_z;
        end
        % -----------------------------------------------------------------
        for Min6 = ( i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main)...
                    +i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                    +i_1(1,i_main) + 11 ):...
                   ( i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main)...
                    +i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                    +i_2(1,i_main) + i_1(1,i_main) + 11 )
            rho_particle(i_main,Min6) = rho_mineral_6_z;
        end
        % -----------------------------------------------------------------
        for Min7 = ( i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main)...
                    +i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                    +i_2(1,i_main) + i_1(1,i_main) + 12 ):...
                   ( i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main)...
                    +i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                    +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 12 )
            rho_particle(i_main,Min7) = rho_mineral_7_z;
        end
        % -----------------------------------------------------------------
        for Min8 = ( i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main)...
                    +i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                    +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 13 ):...
                   ( i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main)...
                    +i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                    +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 13 )
            rho_particle(i_main,Min8) = rho_mineral_8_z;
        end
        % -----------------------------------------------------------------
        for Min9 = ( i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main)...
                    +i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                    +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 14 ):...
                   ( i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main)...
                    +i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main)...
                    +i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 14 )
            rho_particle(i_main,Min9) = rho_mineral_9_z;
        end
        % -----------------------------------------------------------------
        for Min10 = ( i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main)...
                     +i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main)...
                     +i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 15 ):...
                   ( i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main)...
                    +i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main)...
                    +i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                    +i_1(1,i_main) + 15 )
            rho_particle(i_main,Min10) = rho_mineral_10_z;
        end
        % -----------------------------------------------------------------
        for Min11 = ( i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main)...
                     +i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main)...
                     +i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                     +i_1(1,i_main) + 16 ):...
                   ( i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main)...
                    +i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main)...
                    +i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                    +i_2(1,i_main) + i_1(1,i_main) + 16 )
            rho_particle(i_main,Min11) = rho_mineral_11_z;
        end
        % -----------------------------------------------------------------
        for Min12 = ( i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main)...
                     +i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main)...
                     +i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                     +i_2(1,i_main) + i_1(1,i_main) + 17 ):...
                   ( i_18(1,i_main) + i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main)...
                    +i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main)...
                    +i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                    +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 17 )
            rho_particle(i_main,Min12) = rho_mineral_12_z;
        end
        % -----------------------------------------------------------------
        for Min13 = ( i_18(1,i_main) + i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main)...
                    +i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main)...
                    +i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                    +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 18 ):...
                   ( i_19(1,i_main) + i_18(1,i_main) + i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main)...
                    +i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main)...
                    +i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                    +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 18 )
            rho_particle(i_main,Min13) = rho_mineral_13_z;
        end
        % -----------------------------------------------------------------
        for Min14 = ( i_19(1,i_main) + i_18(1,i_main) + i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main)...
                    +i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main)...
                    +i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                    +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main)+ 19 ):...
                   ( i_20(1,i_main) + i_19(1,i_main) + i_18(1,i_main) + i_17(1,i_main) + i_16(1,i_main)...
                    +i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main)...
                    +i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main)...
                    +i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 19 )
            rho_particle(i_main,Min14) = rho_mineral_14_z;
        end
    else
        % -----------------------------------------------------------------
        for Min  = ( i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                    +i_1(1,i_main) + 6 ):...
                   ( i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                    +i_2(1,i_main) + i_1(1,i_main) + 6 )
            rho_particle(i_main,Min) = rho_mineral_z;
        end
        % -----------------------------------------------------------------
    end

end
% -------------------------------------------------------------------------
for i=1:500
    for j=1:size(rho_particle,2)
        if rho_particle(i,j)==0
        rho_particle(i,j) = rho_particle(i,1);
        end 
    end
end
% -------------------------------------------------------------------------
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% -------------------------------------------------------------------------
% --- START OF INITIAL LOOP -----------------------------------------------
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- calculate settling velocity for 1000 particles at each depth --------
% -------------------------------------------------------------------------
num            = 0;
h              = 1;
u              = 1;
w              = 1;
% -------------------------------------------------------------------------
p_aggave(1,1)  = 0.1;
% -------------------------------------------------------------------------
n_num         = 500;                                % number of depth intervals
max_depth_num = 500;                                % max ocean depth [m]
z_num         = linspace(1,max_depth_num,n_num);    % ocean depth grid [m]
% -------------------------------------------------------------------------
num_1            = ones(1,500);
num_1(1,1:3)     = 20;
num_1(1,3:10)    = 20;
num_1(1,10:20)   = 20;
num_1(1,20:30)   = 22;
num_1(1,30:50)   = 22;
num_1(1,50:100)  = 22;
num_1(1,100:160) = 22;
num_1(1,160:200) = 22;
num_1(1,200:240) = 22;
num_1(1,240:280) = 22;
num_1(1,280:320) = 22;
num_1(1,320:360) = 22;
num_1(1,360:500) = 22;
num = interp1(z_num,num_1,z);
% -------------------------------------------------------------------------
for i = 1:n_depth
    n = round(p_aggave(1,h).*n_depth);
    % ---------------------------------------------------------------------
    h = h + 1;
    % ---------------------------------------------------------------------
    if n == 0
        num_particle = 3;
        beta         = ones(1,n_depth);
    else
        num_particle = 1:1:num(i);
    end
    % ---------------------------------------------------------------------
    for jj = 1:n
        power_beta = num_particle(randi(numel(num_particle)));
        beta_init = 10^round(log10(2^power_beta))*rand;
        beta(1,jj) = round(beta_init);
        if beta(1,jj) < 1
            beta(1,jj) = 2;
        end
        
        for kk = (n+1):n_depth
            beta(1,kk) = 1;
        end
    end
    % ---------------------------------------------------------------------
    Beta_f(w,:) = beta;
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    w = w+1;
    % ---------------------------------------------------------------------
    for m = 1:n_depth
        % -----------------------------------------------------------------
        count_org1  = 0;
        count_org2  = 0;
        count_org3  = 0;
        count_org4  = 0;
        count_org5  = 0;
        % -----------------------------------------------------------------
        if Beta_f(i,m) > 15
            count_particle = 15;
        else 
            count_particle = Beta_f(i,m);
        end
        % -----------------------------------------------------------------
        for mm = 1:count_particle
            rho_particle_1 = rho_particle(i,:);
            rho_random = rho_particle_1(randi(numel(rho_particle_1)));
            % -------------------------------------------------------------
            if rho_random == rho_pico_z
                volume_init(1,mm) = m_pico_z(1,i)./rho_pico_z;
                rho_portion(1,mm) = (1/count_particle).*rho_pico_z;
                count_org1        = count_org1 + 1;
            end
            % -------------------------------------------------------------
            if rho_random == rho_cocc_z || rho_random == rho_arag_z
                volume_init(1,mm) = m_cocc_tot_z(1,i)./rho_cocc_z;
                rho_portion(1,mm) = (1/count_particle).*rho_cocc_z;
                count_org2        = count_org2 + 1;
            end
            % -------------------------------------------------------------
            if rho_random == rho_diat_z
                volume_init(1,mm) = m_diat_tot_z(1,i)./rho_diat_z;
                rho_portion(1,mm) = (1/count_particle).*rho_diat_z;
                count_org4        = count_org4 + 1;
            end
            % -------------------------------------------------------------
            if rho_random == rho_dust_z
                volume_init(1,mm) = m_dust_tot_z(1,i)./rho_dust_z;
                rho_portion(1,mm) = (1/count_particle).*rho_dust_z;
            end
            % -------------------------------------------------------------
            if rho_random == rho_dino_z
                volume_init(1,mm) = m_dino_z(1,i)./rho_dino_z;
                rho_portion(1,mm) = (1/count_particle).*rho_dino_z;
                count_org5        = count_org5 + 1;
            end
             % ------------------------------------------------------------
            if strcmp(r_grain,'PSD')
                % ---------------------------------------------------------
                if rho_random == rho_mineral_1_z && rand < K_mineral_1
                    volume_init(1,mm) = m_mineral_tot_1_z(1,i)./rho_mineral_1_z;
                    rho_portion(1,mm) = (1/count_particle).*rho_mineral_1_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_2_z && rand < K_mineral_2
                     volume_init(1,mm) = m_mineral_tot_2_z(1,i)./rho_mineral_2_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_2_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_3_z && rand < K_mineral_3
                     volume_init(1,mm) = m_mineral_tot_3_z(1,i)./rho_mineral_3_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_3_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_4_z && rand < K_mineral_4
                     volume_init(1,mm) = m_mineral_tot_4_z(1,i)./rho_mineral_4_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_4_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_5_z && rand < K_mineral_5
                     volume_init(1,mm) = m_mineral_tot_5_z(1,i)./rho_mineral_5_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_5_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_6_z && rand < K_mineral_6
                     volume_init(1,mm) = m_mineral_tot_6_z(1,i)./rho_mineral_6_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_6_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_7_z && rand < K_mineral_7
                     volume_init(1,mm) = m_mineral_tot_7_z(1,i)./rho_mineral_7_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_7_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_8_z && rand < K_mineral_8
                     volume_init(1,mm) = m_mineral_tot_8_z(1,i)./rho_mineral_8_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_8_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_9_z && rand < K_mineral_9
                     volume_init(1,mm) = m_mineral_tot_9_z(1,i)./rho_mineral_9_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_9_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_10_z && rand < K_mineral_10
                     volume_init(1,mm) = m_mineral_tot_10_z(1,i)./rho_mineral_10_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_10_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_11_z && rand < K_mineral_11
                     volume_init(1,mm) = m_mineral_tot_11_z(1,i)./rho_mineral_11_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_11_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_12_z && rand < K_mineral_12
                     volume_init(1,mm) = m_mineral_tot_12_z(1,i)./rho_mineral_12_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_12_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_13_z && rand < K_mineral_13
                     volume_init(1,mm) = m_mineral_tot_13_z(1,i)./rho_mineral_13_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_13_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_14_z && rand < K_mineral_14
                     volume_init(1,mm) = m_mineral_tot_14_z(1,i)./rho_mineral_14_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_14_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_mineral_14_z
                     volume_init(1,mm) = m_mineral_tot_max_z(1,i)./rho_mineral_14_z;
                     rho_portion(1,mm) = (1/count_particle).*rho_mineral_14_z;
                end
            % -------------------------------------------------------------
            else
            % ---------------------------------------------------------
                if rho_random == rho_mineral_z && rand < K_mineral
                    volume_init(1,mm) = m_mineral_tot_z(1,i)./rho_mineral_z;
                    rho_portion(1,mm) = (1/count_particle).*rho_mineral_z;
                end
            end
            % -------------------------------------------------------------
            volume                  = round(Beta_f(i,m)./count_particle)*sum(volume_init);
            rho_solidaggregate(i,m) = sum(rho_portion);
            % -------------------------------------------------------------  
            % -------------------------------------------------------------
            sum_particles       = count_particle;                          % total number of particles in each aggregate
            count_org_sum       = count_org1 ...
                                 +count_org2 ...
                                 +count_org3 ...
                                 +count_org4 ...
                                 +count_org5;
            C_org_wf(i,m)       = count_org_sum./sum_particles;            % fraction of organic matter particles in aggregate
            C_org_only(i,m)     = (count_org1 + (count_org2*f_org_cocc_z) + (count_org4*f_org_diat) + count_org5)./sum_particles;
            % -------------------------------------------------------------
            volume_particle     = volume./Beta_f(i,m);                     % volume of each particle in aggregate
            r_p(i,m)            = (3*volume_particle./(4*pi)).^0.33333;    % average radius of particles in aggregate [cm]
            r_a(i,m)            = r_p(i,m).* ...
                                  1.e-2* ...
                                  (Beta_f(i,m).^(1./fractal_dim));         % aggregate radius [m]
            % -------------------------------------------------------------
            porosity(i,m)       = 1-Beta_f(i,m).* ...
                                  (((r_p(i,m).*1.e-2)./r_a(i,m)).^3);
            if porosity(i,m) < 0.85
                porosity(i,m)   = 0.85;
            end                                                            % porosity of aggregates
            % -------------------------------------------------------------
            rho_aggregate(i,m)  = (1-porosity(i,m)).* ...
                                  rho_solidaggregate(i,m) ...
                                  +rho_seawater(1,i);                      % density of aggregates [g/cm3]
            % -------------------------------------------------------------
            velocity(i,m)       = ( 0.2222.* ...
                                    (r_a(i,m).^2).* ...
                                    (rho_aggregate(i,m)-rho_seawater(1,i)).* ...
                                    9.81.* ...
                                    7.46.* ...
                                    1.e9 )./ ...
                                  (rho_seawater(1,i).*nho(1,i));           % settling velocity of aggregates [m/day]
            % -------------------------------------------------------------
            Beta(i,m)           = Beta_f(i,m);
            velocity_depth(i,m) = velocity(i,m);
            r_aa(i,m)           = r_a(i,m);
        
        end

    end
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    cc = 1;
    % ---------------------------------------------------------------------
    for t = 1:n_depth
        for k = 1:n_depth
            if k > t
                qq          = min(r_aa(i,t),r_aa(i,k))./max(r_aa(i,t),r_aa(i,k));
                kernel      = ((8.*k_boltz.*(T_seawater(1,i)+273.15))./(6.*mu)).*(((r_aa(i,t)+r_aa(i,k)).^2)./(r_aa(i,t).*r_aa(i,k))+ ...
                              9.8.*(qq^2./(1+2.*qq^2)).*sqrt((epsilon/nho(1,i))).*((r_aa(i,t)+r_aa(i,k)).^3))+ ...
                              0.5.*pi.*min(r_aa(i,t),r_aa(i,k)).^2.*abs(velocity_depth(i,t)-velocity_depth(i,k));
                p_agg(t,cc) = kernel./kernel_ref.*F_NPP_z(i)./F_NPP;
                cc          = cc+1;
            end
        end
    end
    % ---------------------------------------------------------------------
    p_aggave(1,h) = mean(p_agg,'all');
    if p_aggave(1,h) >= 1
        p_aggave(1,h) = 1;
    elseif p_aggave(1,h) < 0.1
        p_aggave(1,h) = 0.1;
    end
    % ---------------------------------------------------------------------
    velocity_average    = mean(velocity_depth(i,:));
    velocity_total(1,u) = velocity_average;
    % ---------------------------------------------------------------------
    u = u+1;
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- ZOOPLANKTON ---------------------------------------------------------
% -------------------------------------------------------------------------
% --- initialize encounter rates ------------------------------------------
% -------------------------------------------------------------------------
F_encounter = (F_NPP./F_NPP_ref).*(z.^-0.2);        % encounter rate b/t aggregates and zooplankton
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- set up fecal pellet size matrix -------------------------------------
% -------------------------------------------------------------------------
fecal1_low = 50; fecal1_high = 100;
fecal2_low = 50; fecal2_high = 100;
fecal3_low = 10; fecal3_high  = 100;
fecal_size = ones(1,n_depth);
% -------------------------------------------------------------------------
alpha_fragm = 1*(z.^-0.65);             % change in fragmentation efficiency with depth
% -------------------------------------------------------------------------
porosity_fecal1_low = 0.5; porosity_fecal1_high = 0.6;
porosity_fecal2_low = 0.5; porosity_fecal2_high = 0.6;
porosity_fecal3_low = 0.5; porosity_fecal3_high = 0.6;
% -------------------------------------------------------------------------
u1 = 1;
% -------------------------------------------------------------------------
for i1 = 1:n_depth
    for m1 = 1:n_depth
        if rand < F_encounter(1,i1)
            Rbreak = 1.e3*atan((r_aa(i1,m1)*1.e6)/1.e4)*alpha_fragm(1,i1);
            if rand < Rbreak
                for i2 = 2:24
                    p_frag(1,i2-1) = 0.91*i2^(-1.56);
                    p_frag_total   = sum(p_frag);
                    if rand < p_frag_total
                        frag        = i2;
                        % -------------------------------------------------
                        r_aa(i1,m1) = r_aa(i1,m1)./frag;                   % fecal pellet radius
                        % -------------------------------------------------
                        velocity(i1,m1) = ( 0.2222.* ...
                                            (r_aa(i1,m1).^2).* ...
                                            (rho_aggregate(i1,m1)-rho_seawater(1,i1)).* ...
                                            9.81.* ...
                                            7.46.* ...
                                            1.e9 )./ ...
                                         (rho_seawater(1,i1).*nho(1,i1));  % settling velocity of fecal pellets [m/day]
                        % -------------------------------------------------
                        velocity_depth(i1,m1) = velocity(i1,m1);
                        Beta(i1,m1)           = round(Beta_f(i1,m1)./frag);
                        if Beta(i1,m1) < 1
                            Beta(i1,m1) = 1;
                        end
                        % -------------------------------------------------
                    end
                end
            elseif rand < C_org_wf(i1,m1)
                % ---------------------------------------------------------
                if i1 <= 40
                    fecal_size = unifrnd(fecal1_low,fecal1_high);
                % ---------------------------------------------------------
                    porosity(i1,m1)      = unifrnd(porosity_fecal1_low,...
                                   porosity_fecal1_high);                  % fecal pellet porosity
                % ---------------------------------------------------------
                elseif i1>40 && i1<=100
                    fecal_size = unifrnd(fecal2_low,fecal2_high);
                % ---------------------------------------------------------
                    porosity(i1,m1)      = unifrnd(porosity_fecal2_low,...
                                   porosity_fecal2_high);                  % fecal pellet porosity
                % ---------------------------------------------------------
                elseif i1>100
                    fecal_size = unifrnd(fecal3_low,fecal3_high);
                % ---------------------------------------------------------
                    porosity(i1,m1)      = unifrnd(porosity_fecal3_low,...
                                   porosity_fecal3_high);                  % fecal pellet porosity
                % ---------------------------------------------------------
                end
                % ---------------------------------------------------------
                r_aa(i1,m1)          = fecal_size.*1.e-6;                  % fecal pellet radius [um]
                % ---------------------------------------------------------
                rho_aggregate(i1,m1) = (1-porosity(i1,m1)).* ...
                                        rho_solidaggregate(i1,m1)+ ...
                                        rho_seawater(1,i1);                % fecal pellet density [g/cm3]
                % ---------------------------------------------------------
                velocity(i1,m1)      = ( 0.2222.* ...
                                         (r_aa(i1,m1).^2).* ...
                                         (rho_aggregate(i1,m1)-rho_seawater(1,i1)).* ...
                                         9.81.* ...
                                         7.46.* ...
                                         1.e9 )./ ...
                                       (rho_seawater(1,i1).*nho(1,i1));    % settling velocity of fecal pellets [m/day]
                % ---------------------------------------------------------
                velocity_depth(i1,m1) = velocity(i1,m1);
            end
        end
    end
    % ---------------------------------------------------------------------
    velocity_average = mean(velocity_depth(i1,:));
    velocity_total(1,u1) = velocity_average;
    % ---------------------------------------------------------------------
    u1 = u1+1;
    % ---------------------------------------------------------------------
end
%
% -------------------------------------------------------------------------
% --- END OF INITIAL LOOP -------------------------------------------------
% -------------------------------------------------------------------------

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %
% -------------------------------------------------------------------------
% --- DOWNSTREAM CALCULATIONS ---------------------------------------------
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- high-res grid for organic matter and silicate calculations ----------
% -------------------------------------------------------------------------
n_hi           = 4.e5;                              % number of depth intervals in high-res grid
z_hi           = linspace(1,max_depth,n_hi);        % high-res depth grid [m]
dz_hi          = max_depth/n_hi;                    % high-res depth interval [m]
velocity_hi    = interp1(z,velocity_total,z_hi);    % interpolated high-res settling velocity [m/d]
K_temperature_hi = interp1(z,K_temperature,z_hi);   % interpolated high-res temperature dependency factor
%
% -------------------------------------------------------------------------
% --- organic matter concentrations and remineralization rates ------------
% -------------------------------------------------------------------------
org_age_d      = cumsum(dz./(velocity_total));                                           % organic matter age [d]
org_top        = (F_NPP/(mean(velocity_total).*365))/12.01;                   % organic matter concentration at top of domain [mol/m3]
org_age_yr     = org_age_d./365;                                              % organic matter age [yr]
    % ---------------------------------------------------------------------
org_age_interp = interp1(z,org_age_yr,z_hi);                                  % interpolaged organic matter age [yr]
k_org          = 10.^(-a_powerlaw*log10(org_age_interp) - b_powerlaw);                  % organic matter reactivity [/yr]
G              = org_top ...
                 .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_org.*dz_hi)); % organic matter concentration at depth [mol/m3] 
% -------------------------------------------------------------------------
Rate_org       = k_org.*G;                                                    % organic matter remineralization rate [mol/m3*yr]
Rate_org_int   = cumsum(Rate_org.*dz_hi);                                     % depth-integrated organic matter remineralization [mol/m2*yr]
BE_org         = G(:)./G(1);                                                  % organic matter burial efficiency at depth
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- input parameters for carbonate system and oxygen --------------------
% -------------------------------------------------------------------------
pO2_mod          = 0.2095;          % modern atmospheric pO2 [atm]
PAL              = 1.0;             % multiplier for model pO2 [times present atmospheric level, PAL]
pO2              = PAL*pO2_mod;     % model atmospheric pO2 [atm]
kO2              = 2.0;             % half-saturation constant for aerobic respiration [uM]
kFe              = 100;             % rate constant for iron oxidation [umol-1 y-1]
Iron             = 0.0002.*ones(1,n_depth); 
Calcium_activity = 0.6;             % activity coefficient for Ca
Calcium          = 10000;           % Ca concentration
k_calcite        = 0.03;            % rate constant for calcite precipitation [uM-1 year-1]
Ksp_ca           = 400000;          % solubility product for calcite [uM^2]
%
% -------------------------------------------------------------------------
% --- saturation calculations ---------------------------------------------
% -------------------------------------------------------------------------
SST_low          = SST;
SST_high         = 4.0;
S                = 35.0;
% -------------------------------------------------------------------------
% constants from Garcia + Gordon [1992], combined fit
% -------------------------------------------------------------------------
a0        = 5.80818;
a1        = 3.20684;
a2        = 4.11890;
a3        = 4.93845;
a4        = 1.01567;
a5        = 1.41575;
b0        = -7.01211*1e-3;
b1        = -7.25958*1e-3;
b2        = -7.93335*1e-3;
b3        = -5.54491*1e-3;
c0        = -1.32412*1e-7;
% -------------------------------------------------------------------------
% scaled temperatures [with T in degC]
% -------------------------------------------------------------------------
Ts_low     = log( (298.15-SST_low)./(273.15+SST_low) );
Ts_high    = log( (298.15-SST_high)./(273.15+SST_high) );
% -------------------------------------------------------------------------
lnC_low    = a0 + a1.*Ts_low + a2.*Ts_low.^2 + a3.*Ts_low.^2 + a3.*Ts_low.^3 + a4.*Ts_low.^4 + a5.*Ts_low.^5 ...
            +S.*(b0 + b1.*Ts_low + b2.*Ts_low.^2 + b3.*Ts_low.^3) ...
            +c0*S^2;
lnC_high   = a0 + a1.*Ts_high + a2.*Ts_high.^2 + a3.*Ts_high.^2 + a3.*Ts_high.^3 + a4.*Ts_high.^4 + a5.*Ts_high.^5 ...
            +S.*(b0 + b1.*Ts_high + b2.*Ts_high.^2 + b3.*Ts_high.^3) ...
            +c0*S^2;
% -------------------------------------------------------------------------
C_O2_low   = exp(lnC_low);
C_O2_high  = exp(lnC_high);
% -------------------------------------------------------------------------
% corrected Henry's Law constants {T,S} [mol L-1 atm-1]
% -------------------------------------------------------------------------
Kh_o2_low   = (C_O2_low/0.2095)*1.e-6;
Kh_o2_high  = (C_O2_high/0.2095)*1.e-6;
% -------------------------------------------------------------------------
% initial (sea surface) concentrations [uM]
% -------------------------------------------------------------------------
O2init      = (Kh_o2_low*pO2)*1.e6;    % low latitude dissolved oxygen
O2init_high = (Kh_o2_high*pO2)*1.e6;   % high latitude dissolved oxygen
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- eddy diffusivity profile --------------------------------------------
% -------------------------------------------------------------------------
Kz_water    = ones(1,n_depth);         % intialize array
Mix_depth   = 200;                     % temperature mixing depth
% -------------------------------------------------------------------------
Kz_water(1,1:Depth_scale/(max_depth/n_depth)) = ((Kz_thermo - Kz_surface)./(Depth_scale)).*(z(1,1:Depth_scale/(max_depth/n_depth))-1) + Kz_surface;
Kz_water(1,(Depth_scale/(max_depth/n_depth)+1):n_depth) =  ((Kz_bottom - Kz_thermo)./(max_depth - Depth_scale)).*...
          (z(1,(Depth_scale/(max_depth/n_depth)+1:n_depth))-z(1,Depth_scale/(max_depth/n_depth)+1)) + Kz_thermo;
Kz_water = Kz_water.*60*60*24*365;     % converted eddy diffusivity [m2 yr-1]
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- depth profile for advective (overturning) fluxes --------------------
% -------------------------------------------------------------------------
K_TH       = 15*1E6;                   % thermohaline overturning [m3/s] (~20 Sv)
A_Downwell = 6*1E13;                   % area of downwelling [m2]
Scaled_TH  = K_TH./A_Downwell;         % scaled overturning [m/s]
Oxygen_TH  = ones(1,n_depth);          % intialize array
TH_depth   = zeros(1,n_depth);         % intialize array
% -------------------------------------------------------------------------
for ii = 1:n_depth
    if z(1,ii) < Depth_scale
        TH_depth(1,ii) = 0;
    else 
        depth_TH (1,ii) = z(1,ii);
        TH_depth(1,ii) = ((2.*Scaled_TH)./(z(1,n_depth)-depth_TH (1,1))).*(z(1,ii)-depth_TH (1,1));
    end
end 
% -------------------------------------------------------------------------
TH_depth = TH_depth.*60*60*24*365;     % converted overturning velocity [m yr-1]; 
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- scaled depth profile for advective (overturning) fluxes -------------
% -------------------------------------------------------------------------
for ii = 1:n_depth 
    if z(1,ii) < Depth_scale
        Oxygen_TH(1,ii) = 0;
    else 
        depth_TH (1,ii) = z(1,ii);
        Oxygen_TH(1,ii) = ((2.*DICinit.*1.5)./(z(1,n_depth)-depth_TH (1,1))).*(z(1,ii)-depth_TH (1,1));     
    end  
end 
% -------------------------------------------------------------------------
Adv_TH  = 1E-10.*TH_depth.*((2.*O2init_high)./(z(1,n_depth)-depth_TH (1,1)));  % scaled advective flux of oxygen [umol/l/year]
Adv_TH1 = K_Adv_TH1.*TH_depth.*((2.*DICinit)./(z(1,n_depth)-depth_TH (1,1)));      % scaled advective flux of DIC [umol/l/year]
Adv_TH22 = 0.01.*TH_depth.*((2.*TAinit)./(z(1,n_depth)-depth_TH (1,1)));      % scaled advective flux of alkalinity [umol/l/year]
% -------------------------------------------------------------------------
if strcmp(Calibration_site,'Pacific')

Adv_TH2 = K1_Adv_TH2.*TH_depth.*((2.*TAinit)./(z(1,n_depth)-depth_TH (1,1)));      % scaled advective flux of alkalinity [umol/l/year]

else 
    
for i = 1:500       
     if (100 < i) && (i < 200)
          Adv_TH2(1,i) = K1_Adv_TH2;
     else
          Adv_TH2(1,i) = K2_Adv_TH2;
     end
end 

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% --- initialize carbonate saturation states ------------------------------
CO3_top          = csys_CO3(TAinit,DICinit);                        % carbonate ion based on DIC and ALK top boundary [uM]
CO3_1            = CO3_top*zeros(1,n_depth);                        % initialize carbonate ion profile
pH_store         = ones(1,n_depth)*csys_pH(TAinit,DICinit);         % initialize pH profile
sigma_carb       = (Calcium_activity.*Calcium.* CO3_1)./Ksp_ca - 1; % intial water column saturation state 
sigma_carb_store = sigma_carb;                                      % store saturation for iteration
CO3_store        = CO3_1;                                           % store carbonate ion for iteration
% -------------------------------------------------------------------------
R_diss_mineral   = 0.0001*ones(1,n_hi);                             % initialize silicate dissolution rate
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
if strcmp(r_grain,'PSD')
    Flux_mineral_z_1_hi  = Flux_mineral_1 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_1.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_1     = interp1(z_hi,Flux_mineral_z_1_hi,z);   
    
    Flux_mineral_z_2_hi  = Flux_mineral_2 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_2.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_2     = interp1(z_hi,Flux_mineral_z_2_hi,z);
    
    Flux_mineral_z_3_hi  = Flux_mineral_3 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_3.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_3     = interp1(z_hi,Flux_mineral_z_3_hi,z);
    
    Flux_mineral_z_4_hi  = Flux_mineral_4 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_4.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_4     = interp1(z_hi,Flux_mineral_z_4_hi,z);
    
    Flux_mineral_z_5_hi  = Flux_mineral_5 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_5.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_5     = interp1(z_hi,Flux_mineral_z_5_hi,z);
    
    Flux_mineral_z_6_hi  = Flux_mineral_6 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_6.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_6     = interp1(z_hi,Flux_mineral_z_6_hi,z);
    
    Flux_mineral_z_7_hi  = Flux_mineral_7 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_7.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_7     = interp1(z_hi,Flux_mineral_z_7_hi,z);
    
    Flux_mineral_z_8_hi  = Flux_mineral_8 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_8.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_8     = interp1(z_hi,Flux_mineral_z_8_hi,z);
    
    Flux_mineral_z_9_hi  = Flux_mineral_9 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_9.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_9     = interp1(z_hi,Flux_mineral_z_9_hi,z);
    
    Flux_mineral_z_10_hi = Flux_mineral_10 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_10.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_10     = interp1(z_hi,Flux_mineral_z_10_hi,z);
    
    Flux_mineral_z_11_hi = Flux_mineral_11 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_11.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_11     = interp1(z_hi,Flux_mineral_z_11_hi,z);
    
    Flux_mineral_z_12_hi = Flux_mineral_12 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_12.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_12     = interp1(z_hi,Flux_mineral_z_12_hi,z);
    
    Flux_mineral_z_13_hi = Flux_mineral_13 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_13.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_13     = interp1(z_hi,Flux_mineral_z_13_hi,z);
    
    Flux_mineral_z_14_hi = Flux_mineral_14 ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_14.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z_14     = interp1(z_hi,Flux_mineral_z_14_hi,z);
    
    Flux_mineral_tot_z =  Flux_mineral_z_1  + Flux_mineral_z_2  + Flux_mineral_z_3 ...
        + Flux_mineral_z_4  + Flux_mineral_z_5  + Flux_mineral_z_6 ...
        + Flux_mineral_z_7  + Flux_mineral_z_8  + Flux_mineral_z_9 ...
        + Flux_mineral_z_10 + Flux_mineral_z_11 + Flux_mineral_z_12 ...
        + Flux_mineral_z_13 + Flux_mineral_z_14;
else
    Flux_mineral_z_hi    = Flux_mineral ...
        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
    Flux_mineral_z     = interp1(z_hi,Flux_mineral_z_hi,z);
end
%
% -------------------------------------------------------------------------
Flux_total_z       =  Flux_pico_z ...
                    + Flux_org_cocc_z ...
                    + Flux_carb_cocc_z ...
                    + Flux_org_arag_z ...
                    + Flux_carb_arag_z ...
                    + Flux_org_diat_z ...
                    + Flux_bSi_diat_z ...
                    + Flux_dino_z ...
                    + Flux_dust_z ...
                    + Flux_mineral_tot_z;                      % total flux
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% --- calculate flux ratios for aggregation probability -------------------
% -------------------------------------------------------------------------
R_pico_z        = Flux_pico_z./Flux_total_z;
R_dino_z        = Flux_dino_z./Flux_total_z;
% -------------------------------------------------------------------------
R_org_cocc_z    = Flux_org_cocc_z./Flux_total_z;
R_carb_cocc_z   = Flux_carb_cocc_z./Flux_total_z;
R_cocc_z        = R_org_cocc_z+R_carb_cocc_z;
% -------------------------------------------------------------------------
R_org_arag_z    = Flux_org_arag_z./Flux_total_z;
R_carb_arag_z   = Flux_carb_arag_z./Flux_total_z;
R_arag_z        = R_org_arag_z+R_carb_arag_z;
% -------------------------------------------------------------------------
R_org_diat_z    = Flux_org_diat_z./Flux_total_z;
R_bSi_diat_z    = Flux_bSi_diat_z./Flux_total_z;
R_diat_z        = R_org_diat_z+R_bSi_diat_z;
% -------------------------------------------------------------------------
R_dust_z        = Flux_dust_z./Flux_total_z;
% -------------------------------------------------------------------------
if strcmp(r_grain,'PSD')
    R_mineral_z_1   = Flux_mineral_z_1./Flux_total_z;
    R_mineral_z_2   = Flux_mineral_z_2./Flux_total_z;
    R_mineral_z_3   = Flux_mineral_z_3./Flux_total_z;
    R_mineral_z_4   = Flux_mineral_z_4./Flux_total_z;
    R_mineral_z_5   = Flux_mineral_z_5./Flux_total_z;
    R_mineral_z_6   = Flux_mineral_z_6./Flux_total_z;
    R_mineral_z_7   = Flux_mineral_z_7./Flux_total_z;
    R_mineral_z_8   = Flux_mineral_z_8./Flux_total_z;
    R_mineral_z_9   = Flux_mineral_z_9./Flux_total_z;
    R_mineral_z_10  = Flux_mineral_z_10./Flux_total_z;
    R_mineral_z_11  = Flux_mineral_z_11./Flux_total_z;
    R_mineral_z_12  = Flux_mineral_z_12./Flux_total_z;
    R_mineral_z_13  = Flux_mineral_z_13./Flux_total_z;
    R_mineral_z_14  = Flux_mineral_z_14./Flux_total_z;
else
    R_mineral_z     = Flux_mineral_z./Flux_total_z;
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
for i=1:500
    
    i_1(1,i)             = round(R_pico_z(1,i).*4.e2);
    i_2(1,i)             = round(R_cocc_z(1,i).*4.e2);
    i_3(1,i)             = round(R_arag_z(1,i).*4.e2);
    i_4(1,i)             = round(R_diat_z(1,i).*4.e2);
    i_5(1,i)             = round(R_dust_z(1,i).*4.e2);
    i_6(1,i)             = round(R_dino_z(1,i).*4.e2);
    % ---------------------------------------------------------------------
    if strcmp(r_grain,'PSD')
        i_7(1,i)             = round(R_mineral_z_1(1,i).*4.e2);
        i_8(1,i)             = round(R_mineral_z_2(1,i).*4.e2);
        i_9(1,i)             = round(R_mineral_z_3(1,i).*4.e2);
        i_10(1,i)            = round(R_mineral_z_4(1,i).*4.e2);
        i_11(1,i)            = round(R_mineral_z_5(1,i).*4.e2);
        i_12(1,i)            = round(R_mineral_z_6(1,i).*4.e2);
        i_13(1,i)            = round(R_mineral_z_7(1,i).*4.e2);
        i_14(1,i)            = round(R_mineral_z_8(1,i).*4.e2);
        i_15(1,i)            = round(R_mineral_z_9(1,i).*4.e2);
        i_16(1,i)            = round(R_mineral_z_10(1,i).*4.e2);
        i_17(1,i)            = round(R_mineral_z_11(1,i).*4.e2);
        i_18(1,i)            = round(R_mineral_z_12(1,i).*4.e2);
        i_19(1,i)            = round(R_mineral_z_13(1,i).*4.e2);
        i_20(1,i)            = round(R_mineral_z_14(1,i).*4.e2);
    else
        i_7(1,i)             = round(R_mineral_z(1,i).*4.e2);
    end
    % ---------------------------------------------------------------------
end
%
for i_main=1:500
    % ---------------------------------------------------------------------
    for i = 1:i_1(1,i_main)
        rho_particle(i_main,i) = rho_pico_z;
    end
    % ---------------------------------------------------------------------
    for j = ( i_1(1,i_main) + 1 ):...
            ( i_2(1,i_main) + i_1(1,i_main) + 1 )
        rho_particle(i_main,j) = rho_cocc_z;
    end
    % ---------------------------------------------------------------------
    for k = ( i_2(1,i_main) + i_1(1,i_main) + 2 ):...
            ( i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 2 )
        rho_particle(i_main,k) = rho_arag_z;
    end
    % ---------------------------------------------------------------------
    for m = ( i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 3 ):...
            ( i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 3 )
        rho_particle(i_main,m) = rho_diat_z;
    end
    % ---------------------------------------------------------------------
    for q = ( i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 4 ):...
            ( i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 4 )
        rho_particle(i_main,q) = rho_dust_z;
    end
    % ---------------------------------------------------------------------
    for pp = ( i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 5 ):...
            ( i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 5 )
        rho_particle(i_main,pp) = rho_dino_z;
    end
    % ---------------------------------------------------------------------
    if strcmp(r_grain,'PSD')
        for Min1 = ( i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                +i_1(1,i_main) + 6 ):...
                ( i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                +i_2(1,i_main) + i_1(1,i_main) + 6 )
            rho_particle(i_main,Min1) = rho_mineral_1_z;
        end
        % -----------------------------------------------------------------
        for Min2 = ( i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                +i_2(1,i_main) + i_1(1,i_main) + 7 ):...
                ( i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 7 )
            rho_particle(i_main,Min2) = rho_mineral_2_z;
        end
        % -----------------------------------------------------------------
        for Min3 = ( i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 8 ):...
                ( i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 8 )
            rho_particle(i_main,Min3) = rho_mineral_3_z;
        end
        % -----------------------------------------------------------------
        for Min4 = ( i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 9 ):...
                ( i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main)...
                +i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 9)
            rho_particle(i_main,Min4) = rho_mineral_4_z;
        end
        % -----------------------------------------------------------------
        for Min5 = ( i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main)...
                +i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 10):...
                ( i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main)...
                +i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                +i_1(1,i_main) + 10 )
            rho_particle(i_main,Min5) = rho_mineral_5_z;
        end
        % -----------------------------------------------------------------
        for Min6 = ( i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main)...
                +i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                +i_1(1,i_main) + 11 ):...
                ( i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main)...
                +i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                +i_2(1,i_main) + i_1(1,i_main) + 11 )
            rho_particle(i_main,Min6) = rho_mineral_6_z;
        end
        % -----------------------------------------------------------------
        for Min7 = ( i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main)...
                +i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                +i_2(1,i_main) + i_1(1,i_main) + 12 ):...
                ( i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main)...
                +i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 12 )
            rho_particle(i_main,Min7) = rho_mineral_7_z;
        end
        % -----------------------------------------------------------------
        for Min8 = ( i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main)...
                +i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 13 ):...
                ( i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main)...
                +i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 13 )
            rho_particle(i_main,Min8) = rho_mineral_8_z;
        end
        % -----------------------------------------------------------------
        for Min9 = ( i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main)...
                +i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 14 ):...
                ( i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main)...
                +i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main)...
                +i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 14 )
            rho_particle(i_main,Min9) = rho_mineral_9_z;
        end
        % -----------------------------------------------------------------
        for Min10 = ( i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main)...
                +i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main)...
                +i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 15 ):...
                ( i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main)...
                +i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main)...
                +i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                +i_1(1,i_main) + 15 )
            rho_particle(i_main,Min10) = rho_mineral_10_z;
        end
        % -----------------------------------------------------------------
        for Min11 = ( i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main)...
                +i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main)...
                +i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                +i_1(1,i_main) + 16 ):...
                ( i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main)...
                +i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main)...
                +i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                +i_2(1,i_main) + i_1(1,i_main) + 16 )
            rho_particle(i_main,Min11) = rho_mineral_11_z;
        end
        % -----------------------------------------------------------------
        for Min12 = ( i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main)...
                +i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main)...
                +i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                +i_2(1,i_main) + i_1(1,i_main) + 17 ):...
                ( i_18(1,i_main) + i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main)...
                +i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main)...
                +i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 17 )
            rho_particle(i_main,Min12) = rho_mineral_12_z;
        end
        % -----------------------------------------------------------------
        for Min13 = ( i_18(1,i_main) + i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main) + i_14(1,i_main)...
                +i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main) + i_9(1,i_main)...
                +i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main)...
                +i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 18 ):...
                ( i_19(1,i_main) + i_18(1,i_main) + i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main)...
                +i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main)...
                +i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 18 )
            rho_particle(i_main,Min13) = rho_mineral_13_z;
        end
        % -----------------------------------------------------------------
        for Min14 = ( i_19(1,i_main) + i_18(1,i_main) + i_17(1,i_main) + i_16(1,i_main) + i_15(1,i_main)...
                +i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main) + i_10(1,i_main)...
                +i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main)...
                +i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main)+ 19 ):...
                ( i_20(1,i_main) + i_19(1,i_main) + i_18(1,i_main) + i_17(1,i_main) + i_16(1,i_main)...
                +i_15(1,i_main) + i_14(1,i_main) + i_13(1,i_main) + i_12(1,i_main) + i_11(1,i_main)...
                +i_10(1,i_main) + i_9(1,i_main) + i_8(1,i_main) + i_7(1,i_main) + i_6(1,i_main)...
                +i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main) + i_1(1,i_main) + 19 )
            rho_particle(i_main,Min14) = rho_mineral_14_z;
        end
    else
        % -----------------------------------------------------------------
        for Min  = ( i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main) + i_2(1,i_main)...
                +i_1(1,i_main) + 6 ):...
                ( i_7(1,i_main) + i_6(1,i_main) + i_5(1,i_main) + i_4(1,i_main) + i_3(1,i_main)...
                +i_2(1,i_main) + i_1(1,i_main) + 6 )
            rho_particle(i_main,Min) = rho_mineral_z;
        end
        % -----------------------------------------------------------------
    end
    % ---------------------------------------------------------------------
end
%
% -------------------------------------------------------------------------
for i=1:500
    for j=1:size(rho_particle,2)
        if rho_particle(i,j)==0
            rho_particle(i,j) = rho_particle(i,1);
        end
    end
end
% -------------------------------------------------------------------------
%
% #########################################################################
% ### START OF PARTICLE FACTORY ITERATION #################################
% #########################################################################
%
% -------------------------------------------------------------------------
% --- set criteria for iteration convergence (main particle factory) ------
% -------------------------------------------------------------------------
K_converge          = 1.0;
iteration_step      = 1.0;
iteration_tolerance = 0.1;
%
% -------------------------------------------------------------------------
% --- set criteria for iteration convergence (carbonate system) -----------
% -------------------------------------------------------------------------
K_converge1           = 1.0;
iteration_step1       = 1.0;
iteration_tolerance_2 = 0.1;
%
% -------------------------------------------------------------------------
% --- iteration for organic matter and silicate dissolution ---------------
% -------------------------------------------------------------------------
while abs(K_converge) > iteration_tolerance || abs(K_converge1) > iteration_tolerance_2
    % ---------------------------------------------------------------------
    % --- initial mass calculations for materials -------------------------
    % ---------------------------------------------------------------------
    % NOTE: all initial masses in units of g/particle
    m_org_pico_z1  = m_org_pico ...
                    *exp(-cumsum((k_org)./(velocity_hi.*365).*dz_hi));  % picoplankton organic mass at depth [g/particle]
    m_org_pico_z   = interp1(z_hi,m_org_pico_z1,z);
    % ---------------------------------------------------------------------
    m_org_cocc_z1  = m_org_cocc ...
                    *exp(-cumsum((k_org)./(velocity_hi.*365).*dz_hi));  % coccolith organic mass at depth [g/particle]
    m_org_cocc_z   = interp1(z_hi,m_org_cocc_z1,z);
    % ---------------------------------------------------------------------
    m_org_arag_z1  = m_org_arag ...
                    *exp(-cumsum((k_org)./(velocity_hi.*365).*dz_hi));  % aragonite-producing phytoplankton organic mass at depth [g/particle]
    m_org_arag_z   = interp1(z_hi,m_org_arag_z1,z);
    % ---------------------------------------------------------------------
    m_org_diat_z1  = m_org_diat ...
                    *exp(-cumsum((k_org)./(velocity_hi.*365).*dz_hi));  % diatom organic mass at depth [g/particle]
    m_org_diat_z   = interp1(z_hi,m_org_diat_z1,z);
    % ---------------------------------------------------------------------
    m_org_dino_z1  = m_org_dino ...
                    *exp(-cumsum((k_org)./(velocity_hi.*365).*dz_hi));  % dinoflagellate organic mass at depth [g/particle]
    m_org_dino_z   = interp1(z_hi,m_org_dino_z1,z);
    % ---------------------------------------------------------------------
    m_carb_cocc_z  = max(m_carb_cocc-k_diss_calc.*z.*1.e-12,num_small); % coccolith calcite mass at depth [g/particle]
    m_carb_arag_z  = max(m_carb_arag-k_diss_arag.*z*1.e-12,num_small);  % aragonite-producing phytoplankton aragonite mass at depth [g/particle]
    m_bSi_diat_z   = ones(1,n_depth).*m_bSi_diat;                       % diatom biogenic silica mass at depth [g/particle]
    % ---------------------------------------------------------------------
    % NOTE: all size classes [n = 1:14] for silicate mineral in units of g/particle
    % ---------------------------------------------------------------------
    if strcmp(r_grain,'PSD')
        m_mineral_1_z1 = m_mineral_1 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_1.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_1_z  = interp1(z_hi,m_mineral_1_z1,z);
        for i=1:n_depth
            if m_mineral_1_z(1,i) < num_small
                m_mineral_1_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------
        m_mineral_2_z1 = m_mineral_2 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_2.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_2_z  = interp1(z_hi,m_mineral_2_z1,z);
        for i=1:n_depth
            if m_mineral_2_z(1,i) < num_small
                m_mineral_2_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------
        m_mineral_3_z1 = m_mineral_3 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_3.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_3_z  = interp1(z_hi,m_mineral_3_z1,z);
        for i=1:n_depth
            if m_mineral_3_z(1,i) < num_small
                m_mineral_3_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------
        m_mineral_4_z1 = m_mineral_4 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_4.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_4_z  = interp1(z_hi,m_mineral_4_z1,z);
        for i=1:n_depth
            if m_mineral_4_z(1,i) < num_small
                m_mineral_4_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------
        m_mineral_5_z1 = m_mineral_5 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_5.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_5_z  = interp1(z_hi,m_mineral_5_z1,z);
        for i=1:n_depth
            if m_mineral_5_z(1,i) < num_small
                m_mineral_5_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------
        m_mineral_6_z1 = m_mineral_6 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_6.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_6_z  = interp1(z_hi,m_mineral_6_z1,z);
        for i=1:n_depth
            if m_mineral_6_z(1,i) < num_small
                m_mineral_6_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------    
        m_mineral_7_z1 = m_mineral_7 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_7.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_7_z  = interp1(z_hi,m_mineral_7_z1,z);
        for i=1:n_depth
            if m_mineral_7_z(1,i) < num_small
                m_mineral_7_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------   
        m_mineral_8_z1 = m_mineral_8 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_8.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_8_z  = interp1(z_hi,m_mineral_8_z1,z);
        for i=1:n_depth
            if m_mineral_8_z(1,i) < num_small
                m_mineral_8_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------    
        m_mineral_9_z1 = m_mineral_9 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_9.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_9_z  = interp1(z_hi,m_mineral_9_z1,z);
        for i=1:n_depth
            if m_mineral_9_z(1,i) < num_small
                m_mineral_9_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------   
        m_mineral_10_z1 = m_mineral_10 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_10.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_10_z  = interp1(z_hi,m_mineral_10_z1,z);
        for i=1:n_depth
            if m_mineral_10_z(1,i) < num_small
                m_mineral_10_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------  
        m_mineral_11_z1 = m_mineral_11 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_11.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_11_z  = interp1(z_hi,m_mineral_11_z1,z);
        for i=1:n_depth
            if m_mineral_11_z(1,i) < num_small
                m_mineral_11_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------   
        m_mineral_12_z1 = m_mineral_12 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_12.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_12_z  = interp1(z_hi,m_mineral_12_z1,z);
        for i=1:n_depth
            if m_mineral_12_z(1,i) < num_small
                m_mineral_12_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------
        m_mineral_13_z1 = m_mineral_13 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_13.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_13_z  = interp1(z_hi,m_mineral_13_z1,z);
        for i=1:n_depth
            if m_mineral_13_z(1,i) < num_small
                m_mineral_13_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------   
        m_mineral_14_z1 = m_mineral_14 ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_14.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_14_z  = interp1(z_hi,m_mineral_14_z1,z);
        for i=1:n_depth
            if m_mineral_14_z(1,i) < num_small
                m_mineral_14_z(1,i) = num_small;
            end
        end
    else
        % -----------------------------------------------------------------
        m_mineral_z1   = m_mineral ...
                        .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral.*K_temperature_hi.*ones(1,n_hi).*dz_hi));
        m_mineral_z    = interp1(z_hi,m_mineral_z1,z);
        for i=1:n_depth
            if m_mineral_z(1,i) < num_small
                m_mineral_z(1,i) = num_small;
            end
        end
        % -----------------------------------------------------------------
    end
    % ---------------------------------------------------------------------
    %
    % ---------------------------------------------------------------------
    % --- determining grain size of most probable size class --------------
    % ---------------------------------------------------------------------
    if strcmp(r_grain,'PSD')
        max_mineral = max([K_mineral_1  K_mineral_2  K_mineral_3  ...
                           K_mineral_4  K_mineral_5  K_mineral_6 ...
                           K_mineral_7  K_mineral_8  K_mineral_9 ...
                           K_mineral_10 K_mineral_11 K_mineral_12 ...
                           K_mineral_13 K_mineral_14]);
        % -----------------------------------------------------------------
        if K_mineral_1 == max_mineral
            m_mineral_max_z = m_mineral_1_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_2 == max_mineral
            m_mineral_max_z = m_mineral_2_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_3 == max_mineral
            m_mineral_max_z = m_mineral_3_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_4 == max_mineral
            m_mineral_max_z = m_mineral_4_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_5 == max_mineral
            m_mineral_max_z = m_mineral_5_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_6 == max_mineral
            m_mineral_max_z = m_mineral_6_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_7 == max_mineral
            m_mineral_max_z = m_mineral_7_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_8 == max_mineral
            m_mineral_max_z = m_mineral_8_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_9 == max_mineral
            m_mineral_max_z = m_mineral_9_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_10 == max_mineral
            m_mineral_max_z = m_mineral_10_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_11 == max_mineral
            m_mineral_max_z = m_mineral_11_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_12 ==  max_mineral
            m_mineral_max_z = m_mineral_12_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_13 == max_mineral
            m_mineral_max_z = m_mineral_13_z;
        end
        % -----------------------------------------------------------------
        if K_mineral_14 == max_mineral
            m_mineral_max_z = m_mineral_14_z;
        end
        % -----------------------------------------------------------------
    else
        % -----------------------------------------------------------------
        max_mineral     = K_mineral;
        m_mineral_max_z = m_mineral_z;
        % -----------------------------------------------------------------
    end
    %
    % ---------------------------------------------------------------------
    F_NPP_z = F_NPP.*exp(-cumsum((k_org)./(velocity_hi.*365).*dz_hi));
    % ---------------------------------------------------------------------
    m_pico_z            = m_org_pico_z;                  % total picoplankton mass
    m_dino_z            = m_org_dino_z;                  % total dinoflagellate mass
    m_cocc_tot_z        = m_org_cocc_z + m_carb_cocc_z;  % total coccolith mass
    m_arag_tot_z        = m_org_arag_z + m_carb_arag_z;  % total aragonite plankton mass
    m_diat_tot_z        = m_org_diat_z + m_bSi_diat_z;   % total diatom mass
    % ---------------------------------------------------------------------
    m_dust_tot_z        = ones(1,n_depth).*m_dust;           % total dust mass
    % ---------------------------------------------------------------------
    % NOTE: mass of silicate in each size class
    if strcmp(r_grain,'PSD')
        m_mineral_tot_1_z   = m_mineral_1_z;
        m_mineral_tot_2_z   = m_mineral_2_z;
        m_mineral_tot_3_z   = m_mineral_3_z;
        m_mineral_tot_4_z   = m_mineral_4_z;
        m_mineral_tot_5_z   = m_mineral_5_z;
        m_mineral_tot_6_z   = m_mineral_6_z;
        m_mineral_tot_7_z   = m_mineral_7_z;
        m_mineral_tot_8_z   = m_mineral_8_z;
        m_mineral_tot_9_z   = m_mineral_9_z;
        m_mineral_tot_10_z  = m_mineral_10_z;
        m_mineral_tot_11_z  = m_mineral_11_z;
        m_mineral_tot_12_z  = m_mineral_12_z;
        m_mineral_tot_13_z  = m_mineral_13_z;
        m_mineral_tot_14_z  = m_mineral_14_z;
        m_mineral_tot_max_z = m_mineral_max_z;
    else
        m_mineral_tot_z     = m_mineral_z;
        m_mineral_tot_max_z = m_mineral_max_z;
    end
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    % --- initialize volumes and densities for particle materials ---------
    % ---------------------------------------------------------------------
    % NOTE: fractions of each material [org|carb|bSi] in plankton
    f_org_cocc_z     = m_org_cocc_z/(m_org_cocc_z+m_carb_cocc_z);
    f_carb_cocc_z    = m_carb_cocc_z/(m_org_cocc_z+m_carb_cocc_z);
    f_org_arag_z     = m_org_arag_z/(m_org_arag_z+m_carb_arag_z);
    f_carb_arag_z    = m_carb_arag_z/(m_org_arag_z+m_carb_arag_z);
    f_org_diat_z     = m_org_diat_z/(m_org_diat_z+m_bSi_diat_z);
    f_bSi_diat_z     = m_bSi_diat_z/(m_org_diat_z+m_bSi_diat_z);
    % ---------------------------------------------------------------------
    % NOTE: densities of plankton calculated via mass balance
    rho_pico_z       = rho_org;
    rho_dino_z       = rho_dino;
    rho_cocc_z       = (f_org_cocc_z*rho_org)+(f_carb_cocc_z*rho_carb_cocc);
    rho_arag_z       = (f_org_arag_z*rho_org)+(f_carb_arag_z*rho_carb_arag);
    rho_diat_z       = (f_org_diat_z*rho_org)+(f_bSi_diat_z*rho_bSi);
    % ---------------------------------------------------------------------
    rho_dust_z       = rho_dust;
    % ---------------------------------------------------------------------
    % NOTE: densities of all size classes defined by mineral phase
    if strcmp(r_grain,'PSD')
        rho_mineral_1_z  = rho_mineral;
        rho_mineral_2_z  = rho_mineral;
        rho_mineral_3_z  = rho_mineral;
        rho_mineral_4_z  = rho_mineral;
        rho_mineral_5_z  = rho_mineral;
        rho_mineral_6_z  = rho_mineral;
        rho_mineral_7_z  = rho_mineral;
        rho_mineral_8_z  = rho_mineral;
        rho_mineral_9_z  = rho_mineral;
        rho_mineral_10_z = rho_mineral;
        rho_mineral_11_z = rho_mineral;
        rho_mineral_12_z = rho_mineral;
        rho_mineral_13_z = rho_mineral;
        rho_mineral_14_z = rho_mineral;
    else
        rho_mineral_z    = rho_mineral;
    end
    % ---------------------------------------------------------------------
    % NOTE: volumes of particle constituents [cm3]
    v_pico_z        = m_pico_z/rho_pico_z;
    v_cocc_z        = m_cocc_tot_z/rho_cocc_z;
    v_arag_z        = m_arag_tot_z/rho_arag_z;
    v_diat_z        = m_diat_tot_z/rho_diat_z;
    v_dino_z        = m_dino_z/rho_dino_z;
    % ---------------------------------------------------------------------
    v_dust_z        = m_dust_tot_z/rho_dust_z;
    % ---------------------------------------------------------------------
    if strcmp(r_grain,'PSD')
        v_mineral_1_z   = m_mineral_tot_1_z/rho_mineral_1_z;
        v_mineral_2_z   = m_mineral_tot_2_z/rho_mineral_2_z;
        v_mineral_3_z   = m_mineral_tot_3_z/rho_mineral_3_z;
        v_mineral_4_z   = m_mineral_tot_4_z/rho_mineral_4_z;
        v_mineral_5_z   = m_mineral_tot_5_z/rho_mineral_5_z;
        v_mineral_6_z   = m_mineral_tot_6_z/rho_mineral_6_z;
        v_mineral_7_z   = m_mineral_tot_7_z/rho_mineral_7_z;
        v_mineral_8_z   = m_mineral_tot_8_z/rho_mineral_8_z;
        v_mineral_9_z   = m_mineral_tot_9_z/rho_mineral_9_z;
        v_mineral_10_z  = m_mineral_tot_10_z/rho_mineral_10_z;
        v_mineral_11_z  = m_mineral_tot_11_z/rho_mineral_11_z;
        v_mineral_12_z  = m_mineral_tot_12_z/rho_mineral_12_z;
        v_mineral_13_z  = m_mineral_tot_13_z/rho_mineral_13_z;
        v_mineral_14_z  = m_mineral_tot_14_z/rho_mineral_14_z;
        v_mineral_max_z = m_mineral_tot_max_z/rho_mineral_14_z;
    else
        v_mineral_z     = m_mineral_tot_z/rho_mineral_z;
        v_mineral_max_z = m_mineral_tot_max_z/rho_mineral_z;
    end
    %
    % ---------------------------------------------------------------------

    % ---------------------------------------------------------------------
    %
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %
    % ---------------------------------------------------------------------
    % --- START OF INITIAL LOOP -------------------------------------------
    % ---------------------------------------------------------------------
    u = 1;
    % ---------------------------------------------------------------------
    % --- calculate settling velocity for 1000 particles at each depth ----
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    for i = 1:n_depth
        % -----------------------------------------------------------------
        for m = 1:n_depth
            % -------------------------------------------------------------
            count_org1  = 0;
            count_org2  = 0;
            count_org3  = 0;
            count_org4  = 0;
            count_org5  = 0;
            % -------------------------------------------------------------
            if Beta_f(i,m) > 15
                count_particle = 15;
            else 
                count_particle = Beta_f(i,m);
            end
            % -------------------------------------------------------------
            for mm = 1:count_particle
                rho_particle_1 = rho_particle(i,:);
                rho_random = rho_particle_1(randi(numel(rho_particle_1)));
                % ---------------------------------------------------------
                if rho_random == rho_pico_z
                    volume_init(1,mm) = m_pico_z(1,i)./rho_pico_z;
                    rho_portion(1,mm) = (1/count_particle).*rho_pico_z;
                    count_org1        = count_org1 + 1;
                end
                % ---------------------------------------------------------
                if rho_random == rho_cocc_z
                    volume_init(1,mm) = m_cocc_tot_z(1,i)./rho_cocc_z;
                    rho_portion(1,mm) = (1/count_particle).*rho_cocc_z;
                    count_org2        = count_org2 + 1;
                end
                % ---------------------------------------------------------
                if rho_random == rho_diat_z
                    volume_init(1,mm) = m_diat_tot_z(1,i)./rho_diat_z;
                    rho_portion(1,mm) = (1/count_particle).*rho_diat_z;
                    count_org4        = count_org4 + 1;
                end
                % ---------------------------------------------------------
                if rho_random == rho_dust_z
                    volume_init(1,mm) = m_dust_tot_z(1,i)./rho_dust_z;
                    rho_portion(1,mm) = (1/count_particle).*rho_dust_z;
                end
                % ---------------------------------------------------------
                if rho_random == rho_dino_z
                    volume_init(1,mm) = m_dino_z(1,i)./rho_dino_z;
                    rho_portion(1,mm) = (1/count_particle).*rho_dino_z;
                    count_org5        = count_org5 + 1;
                end
                 % --------------------------------------------------------
                if strcmp(r_grain,'PSD')
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_1_z && rand < K_mineral_1
                        volume_init(1,mm) = m_mineral_tot_1_z(1,i)./rho_mineral_1_z;
                        rho_portion(1,mm) = (1/count_particle).*rho_mineral_1_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_2_z && rand < K_mineral_2
                         volume_init(1,mm) = m_mineral_tot_2_z(1,i)./rho_mineral_2_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_2_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_3_z && rand < K_mineral_3
                         volume_init(1,mm) = m_mineral_tot_3_z(1,i)./rho_mineral_3_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_3_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_4_z && rand < K_mineral_4
                         volume_init(1,mm) = m_mineral_tot_4_z(1,i)./rho_mineral_4_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_4_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_5_z && rand < K_mineral_5
                         volume_init(1,mm) = m_mineral_tot_5_z(1,i)./rho_mineral_5_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_5_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_6_z && rand < K_mineral_6
                         volume_init(1,mm) = m_mineral_tot_6_z(1,i)./rho_mineral_6_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_6_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_7_z && rand < K_mineral_7
                         volume_init(1,mm) = m_mineral_tot_7_z(1,i)./rho_mineral_7_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_7_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_8_z && rand < K_mineral_8
                         volume_init(1,mm) = m_mineral_tot_8_z(1,i)./rho_mineral_8_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_8_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_9_z && rand < K_mineral_9
                         volume_init(1,mm) = m_mineral_tot_9_z(1,i)./rho_mineral_9_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_9_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_10_z && rand < K_mineral_10
                         volume_init(1,mm) = m_mineral_tot_10_z(1,i)./rho_mineral_10_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_10_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_11_z && rand < K_mineral_11
                         volume_init(1,mm) = m_mineral_tot_11_z(1,i)./rho_mineral_11_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_11_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_12_z && rand < K_mineral_12
                         volume_init(1,mm) = m_mineral_tot_12_z(1,i)./rho_mineral_12_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_12_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_13_z && rand < K_mineral_13
                         volume_init(1,mm) = m_mineral_tot_13_z(1,i)./rho_mineral_13_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_13_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_14_z && rand < K_mineral_14
                         volume_init(1,mm) = m_mineral_tot_14_z(1,i)./rho_mineral_14_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_14_z;
                    end
                    % -----------------------------------------------------
                    if rho_random == rho_mineral_14_z
                         volume_init(1,mm) = m_mineral_tot_max_z(1,i)./rho_mineral_14_z;
                         rho_portion(1,mm) = (1/count_particle).*rho_mineral_14_z;
                    end
                % ---------------------------------------------------------
                else
                % ---------------------------------------------------------
                    if rho_random == rho_mineral_z && rand < K_mineral
                        volume_init(1,mm) = m_mineral_tot_z(1,i)./rho_mineral_z;
                        rho_portion(1,mm) = (1/count_particle).*rho_mineral_z;
                    end
                end
                % ---------------------------------------------------------
                volume                  = round(Beta_f(i,m)./count_particle)*sum(volume_init);
                rho_solidaggregate(i,m) = sum(rho_portion);
                % ---------------------------------------------------------  
                % ---------------------------------------------------------
                sum_particles       = count_particle;                          % total number of particles in each aggregate
                count_org_sum       = count_org1 ...
                                     +count_org2 ...
                                     +count_org3 ...
                                     +count_org4 ...
                                     +count_org5;
                C_org_wf(i,m)       = count_org_sum./sum_particles;            % fraction of organic matter particles in aggregate
                C_org_only(i,m)     = (count_org1 + (count_org2*f_org_cocc_z) + (count_org4*f_org_diat) + count_org5)./sum_particles;
                % ---------------------------------------------------------
                volume_particle     = volume./Beta_f(i,m);                     % volume of each particle in aggregate
                r_p(i,m)            = (3*volume_particle./(4*pi)).^0.33333;    % average radius of particles in aggregate [cm]
                r_a(i,m)            = r_p(i,m).* ...
                                      1.e-2* ...
                                      (Beta_f(i,m).^(1./fractal_dim));         % aggregate radius [m]
                % ---------------------------------------------------------
                porosity(i,m)       = 1-Beta_f(i,m).* ...
                                      (((r_p(i,m).*1.e-2)./r_a(i,m)).^3);
                if porosity(i,m) < 0.85
                    porosity(i,m)   = 0.85;
                end                                                            % porosity of aggregates
                % ---------------------------------------------------------
                rho_aggregate(i,m)  = (1-porosity(i,m)).* ...
                                      rho_solidaggregate(i,m) ...
                                      +rho_seawater(1,i);                      % density of aggregates [g/cm3]
                % ---------------------------------------------------------
                velocity(i,m)       = ( 0.2222.* ...
                                        (r_a(i,m).^2).* ...
                                        (rho_aggregate(i,m)-rho_seawater(1,i)).* ...
                                        9.81.* ...
                                        7.46.* ...
                                        1.e9 )./ ...
                                      (rho_seawater(1,i).*nho(1,i));           % settling velocity of aggregates [m/day]
                % ---------------------------------------------------------
                Beta(i,m)           = Beta_f(i,m);
                velocity_depth(i,m) = velocity(i,m);
                r_aa(i,m)           = r_a(i,m);

            end

        end
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        velocity_average    = mean(velocity_depth(i,:));
        velocity_total(1,u) = velocity_average;
        % -----------------------------------------------------------------
        u = u+1;
        % -----------------------------------------------------------------
    end
    % ---------------------------------------------------------------------
    %
    % ---------------------------------------------------------------------
    % --- ZOOPLANKTON -----------------------------------------------------
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %
    % ---------------------------------------------------------------------
    % --- set up fecal pellet size matrix ---------------------------------
    % ---------------------------------------------------------------------
    u1 = 1;
    count_zoo_i = 1;
    count_fecal_i = 1;
    count_fecal_m = 1;
    % ---------------------------------------------------------------------
    for i1 = 1:n_depth
        count_zoo_m = 1;
        count_fecal_m = 1;
        for m1 = 1:n_depth
            if rand < F_encounter(1,i1)
                Rbreak = 1.e3*atan((r_aa(i1,m1)*1.e6)/1.e4)*alpha_fragm(1,i1);
                count_encounter_zoo(count_zoo_i,count_zoo_m) = 1;
                count_zoo_m = count_zoo_m + 1;
                if rand < Rbreak
                    for i2 = 2:24
                        p_frag(1,i2-1) = 0.91*i2^(-1.56);
                        p_frag_total   = sum(p_frag);
                        if rand < p_frag_total
                            frag        = i2;
                            % ---------------------------------------------
                            r_aa(i1,m1) = r_aa(i1,m1)./frag;                   % fecal pellet radius
                            % ---------------------------------------------
                            velocity(i1,m1) = ( 0.2222.* ...
                                                (r_aa(i1,m1).^2).* ...
                                                (rho_aggregate(i1,m1)-rho_seawater(1,i1)).* ...
                                                9.81.* ...
                                                7.46.* ...
                                                1.e9 )./ ...
                                             (rho_seawater(1,i1).*nho(1,i1));  % settling velocity of fecal pellets [m/day]
                            % ---------------------------------------------
                            velocity_depth(i1,m1) = velocity(i1,m1);
                            Beta(i1,m1)           = round(Beta_f(i1,m1)./frag);
                            if Beta(i1,m1) < 1
                                Beta(i1,m1) = 1;
                            end
                            % ---------------------------------------------
                        end
                    end
                elseif rand < C_org_wf(i1,m1)
                    % -----------------------------------------------------
                    count_fecal_zoo(count_fecal_i,count_fecal_m) = 1;
                    count_fecal_m = count_fecal_m + 1;
                    % -----------------------------------------------------
                    if i1 <= 20
                        fecal_size = unifrnd(fecal1_low,fecal1_high);
                    % -----------------------------------------------------
                        porosity(i1,m1)   = unifrnd(porosity_fecal1_low,...
                                       porosity_fecal1_high);                   % fecal pellet porosity
                    % -----------------------------------------------------
                    elseif i1>20 && i1<=100
                        fecal_size = unifrnd(fecal2_low,fecal2_high);
                    % -----------------------------------------------------
                        porosity(i1,m1)   = unifrnd(porosity_fecal2_low,...
                                       porosity_fecal2_high);                   % fecal pellet porosity
                    % -----------------------------------------------------
                    elseif i1>100
                        fecal_size = unifrnd(fecal3_low,fecal3_high);
                    % -----------------------------------------------------
                        porosity(i1,m1)   = unifrnd(porosity_fecal3_low,...
                                       porosity_fecal3_high);                   % fecal pellet porosity
                    % -----------------------------------------------------
                    end
                    % -----------------------------------------------------
                    r_aa(i1,m1)          = fecal_size.*1.e-6;                   % fecal pellet radius [um]
                    % -----------------------------------------------------
                    rho_aggregate(i1,m1) = (1-porosity(i1,m1)).* ...
                                            rho_solidaggregate(i1,m1)+ ...
                                            rho_seawater(1,i1);                 % fecal pellet density [g/cm3]
                    % -----------------------------------------------------
                    velocity(i1,m1)      = ( 0.2222.* ...
                                             (r_aa(i1,m1).^2).* ...
                                             (rho_aggregate(i1,m1)-rho_seawater(1,i1)).* ...
                                             9.81.* ...
                                             7.46.* ...
                                             1.e9 )./ ...
                                           (rho_seawater(1,i1).*nho(1,i1));     % settling velocity of fecal pellets [m/day]
                    % -----------------------------------------------------
                    velocity_depth(i1,m1) = velocity(i1,m1);
                end
            end
        end
        % ---------------------------------------------------------------------
        velocity_average = mean(velocity_depth(i1,:));
        velocity_total(1,u1) = velocity_average;
        % -----------------------------------------------------------------
        u1 = u1+1;
        count_zoo_i = count_zoo_i + 1;
        count_fecal_i = count_fecal_i + 1;
        % -----------------------------------------------------------------
    end
    % --------------------------------------------------------------------- 
    Percentage_encounter_zoo = sum(count_encounter_zoo')./n_depth;
    Percentage_fecal_zoo = sum(count_fecal_zoo')./n_depth;
    % ---------------------------------------------------------------------
    
    % #####################################################################
    % ### END OF MAIN LOOP ################################################
    % #####################################################################
    %
    % ---------------------------------------------------------------------
    % --- high-res grid for organic matter and silicate calculations ------
    % ---------------------------------------------------------------------
    n_hi           = 4.e5;                              % number of depth intervals in high-res grid
    z_hi           = linspace(1,max_depth,n_hi);        % high-res depth grid [m]
    dz_hi          = max_depth/n_hi;                    % high-res depth interval [m]
    velocity_hi    = interp1(z,velocity_total,z_hi);    % interpolated high-res settling velocity [m/d]
    %
    % ---------------------------------------------------------------------
    % --- organic matter concentrations and remineralization rates --------
    % ---------------------------------------------------------------------
    org_age_d      = cumsum(dz./velocity_total);                                          % organic matter age [d]
    org_top        = (F_NPP/(mean(velocity_total).*365))/12.01;                  % organic matter concentration at top of domain [mol/m3]
    org_age_yr     = org_age_d./365;                                             % organic matter age [yr]
    % ---------------------------------------------------------------------
    org_age_interp = interp1(z,org_age_yr,z_hi);                                 % interpolaged organic matter age [yr]
    k_org          = 10.^(-a_powerlaw*log10(org_age_interp) - b_powerlaw);                 % organic matter reactivity [/yr]
    G              = org_top ...
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_org.*dz_hi)); % organic matter concentration at depth [mol/m3]
    % ---------------------------------------------------------------------
    Rate_org       = k_org.*G;                                                   % organic matter remineralization rate [mol/m2*yr]
    Rate_org_int   = cumsum(Rate_org.*dz_hi);                                    % depth-integrated organic matter remineralization [mol/m2*yr]
    BE_org         = G(:)./G(1);                                                 % organic matter burial efficiency at depth
    % ---------------------------------------------------------------------
    %
    % ---------------------------------------------------------------------
    % -- step through iteration -------------------------------------------
    % ---------------------------------------------------------------------
    iteration_step                       = iteration_step + 1;
    Velocity_iteration(iteration_step,:) = velocity_total;
    % ---------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % --- evaluate convergence criterion for main particle factory --------
    % ---------------------------------------------------------------------
    K_converge  = ( mean(Velocity_iteration(iteration_step,:)) ...
                   -mean(Velocity_iteration(iteration_step-1,:)) )./...
                  (mean(Velocity_iteration(iteration_step-1,:)));
    % ---------------------------------------------------------------------
    %                                                      
    % ---------------------------------------------------------------------
    % --- solve ODE for dissolved oxygen ----------------------------------
    % ---------------------------------------------------------------------
    Rate_org_interp  = interp1(z_hi,Rate_org,z);                           % interpolated remineralization rate [mol/m3/yr]
    % ---------------------------------------------------------------------
    n_O2             = n_depth;
    nmesh            = 100;
    solinit          = bvpinit(linspace(1,max_depth,nmesh),[1 0]);
    sol              = bvp4c(@ODE_O2,@BC_O2,solinit);
    x                = linspace(1,max_depth,n_O2);
    p                = deval(sol,x);
    Oxygen           = p(1,:);
    % ---------------------------------------------------------------------
    Rate_O2_resp     = Rate_org_interp.*Oxygen./(Oxygen+kO2);  % rate of aerobic respiration [mol/m3/yr]
    % ---------------------------------------------------------------------
    %
    % ---------------------------------------------------------------------
    % --- carbonate precipitation/dissolution  profile --------------------
    % ---------------------------------------------------------------------
    for i=1:n_depth
        if sigma_carb(1,i) > 1
            K_sigma      =   abs(sigma_carb(1,i)-1).^1.7;
            R1_carb(1,i) =   K_sigma.*k_calcite;        % carbonate precipitation [umol/l/yr]
        elseif sigma_carb(1,i) < 1
            K_sigma      =   abs(sigma_carb(1,i)-1).^1.7;
            R1_carb(1,i) = - K_sigma.*k_calcite;        % carbonate dissolution [umol/l/yr]
        else
            R1_carb(1,i) = 0.0;
        end
    end  
    % ---------------------------------------------------------------------
    %
    % ---------------------------------------------------------------------
    % --- solve DIC profile -----------------------------------------------
    % ---------------------------------------------------------------------
    n_DIC            = n_depth;
    nmesh            = 100;
    solinit          = bvpinit(linspace(1,max_depth,nmesh),[1 0]); 
    sol              = bvp4c(@ODE_DIC,@BC_DIC,solinit);
    x                = linspace(1,max_depth,n_DIC);
    DIC              = deval(sol,x);
    C_DIC            = DIC(1,:);
    % ---------------------------------------------------------------------
    %
    % ---------------------------------------------------------------------
    % --- solve alkalinity profile ----------------------------------------
    % ---------------------------------------------------------------------
    R_diss_mineral_TA = interp1(z_hi,R_diss_mineral,z);  % interpolated high-res dissolution profile 
    R_diss_mineral_1   = R_diss_mineral_TA*1.e3;         % mineral dissolution rate [umol/l/y]
    % ---------------------------------------------------------------------
    n_TA              = n_depth;
    nmesh             = 100;
    solinit           = bvpinit(linspace(1,max_depth,nmesh),[1 0]); 
    sol               = bvp4c(@ODE_TA,@BC_TA,solinit);
    x                 = linspace(1,max_depth,n_TA);
    TA_conc           = deval(sol,x);
    ALK               = TA_conc(1,:);
    % ---------------------------------------------------------------------
    %
    % --- pH calculation --------------------------------------------------
    for i=1:n_depth
        pH_1(1,i)    = csys_pH(ALK(1,i),C_DIC(1,i));
        C_H2CO3(1,i) = csys_H2CO3(ALK(1,i),C_DIC(1,i));
        CO3_1(1,i)   = csys_CO3(ALK(1,i),C_DIC(1,i));
    end
    % ---------------------------------------------------------------------
    pH = pH_1;
    % ---------------------------------------------------------------------
    %
    % ---------------------------------------------------------------------
    % --- correcting dissolution rates for pH -----------------------------
    % ---------------------------------------------------------------------
    pH_hi    = interp1(z,pH,z_hi);    % interpolated high-res pH profile
    T_seawater_hi = interp1(z,T_seawater,z_hi); % interpolated high-res temperature profile
    % ---------------------------------------------------------------------
    if strcmp(mineral,'basalt')
        %k_min_pH0 = p1.*(pH_hi.^3) + p2.*(pH_hi.^2) + p3.*pH_hi + p4;
        k_min_pH0 = log10(mineral_mass.*(A_0a*exp(-E_aa./(R.*(T_seawater_hi+C_to_K))).*10.^(-n_a.*pH_hi) + A_0n*exp(-E_an./(R.*(T_seawater_hi+C_to_K))).*10.^(-n_n.*pH_hi) + A_0b*exp(-E_ab./(R.*(T_seawater_hi+C_to_K))).*10.^(-n_b.*pH_hi))); % see Gudbrandsson et al. [2011] | doi:10.1016/j.gca.2011.06.035
    end
    if strcmp(mineral,'olivine')
        %k_min_pH0 = (p1.*pH_hi + p2);
        k_min_pH0 = log10(mineral_mass.*(A_0a*exp(-E_aa./(R.*(T_seawater_hi+C_to_K))).*10.^(-n_a.*pH_hi) + A_0b*exp(-E_ab./(R.*(T_seawater_hi+C_to_K))).*10.^(-n_b.*pH_hi)));   % see Oelkers et al. [2018] | doi:10.1016/j.chemgeo.2018.10.008
    end
    % ---------------------------------------------------------------------
    if strcmp(mineral,'MgO')
        k_min_pH0 = log10(mineral_mass.*(A_0a*exp(-E_aa./(R.*(T_seawater_hi+C_to_K))).*10.^(-n_a.*pH_hi) + A_0n*exp(-E_an./(R.*(T_seawater_hi+C_to_K))).*10.^(-n_n.*pH_hi)));   % see Pokrovsky and Schott. [2004] | doi:10.1016/S0016-7037(03)00238-2
    end
    % ---------------------------------------------------------------------
    if strcmp(mineral,'CaO')
        k_min_pH0 = log10(mineral_mass.*(A_0*exp(-E_a./(R.*(T_seawater_hi+C_to_K))).*10.^(-n_all.*pH_hi))); % see Wang et al. [1998] | doi:10.1080/00986449808912726
    end
    % ---------------------------------------------------------------------
    if 7 < pH(end) && pH(end) < 10 || iteration_step1 == 1 
        iteration_step1 = iteration_step1 + 1;
        pH_store(iteration_step1,:) = pH;
    end
    % ---------------------------------------------------------------------
    %
    % --- update carbonate saturation states ------------------------------
    sigma_carb                          = (Calcium_activity.*Calcium.* CO3_1)./Ksp_ca - 1;  
    sigma_carb_store(iteration_step1,:) = sigma_carb;
    CO3_store(iteration_step1,:)        = CO3_1;
    % ---------------------------------------------------------------------
    %
    % --- update silicate dissolution fluxes ------------------------------
    if strcmp(r_grain,'PSD')
        % -----------------------------------------------------------------
        % --- size class 1 ------------------------------------------------
        % -----------------------------------------------------------------  
        k_sed_min_1 = SSA_1.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                    %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_1 = k_sed_min_1.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));   
        F_mineral_11  =  Flux_mineral_1 / mineral_mass;
        F_mineral_1_z   = F_mineral_11 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_1.*K_temperature_hi.*dz_hi) ...
                        );                 
        mineral_1_conc  = (F_mineral_11./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_1.*K_temperature_hi.*dz_hi) ...
                        );                                                          % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral1 = k_diss_mineral_1.* K_temperature_hi.* mineral_1_conc;     % mineral dissolution rate [mol/m3*y]
        R_diss_int_1     = cumsum(R_diss_mineral1.*dz_hi);                          % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_1     = mineral_1_conc(:)/mineral_1_conc(1);                     % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 2 ------------------------------------------------
        % -----------------------------------------------------------------     
        k_sed_min_2 = SSA_2.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                    %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_2 = k_sed_min_2.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));   
        F_mineral_22  =  Flux_mineral_2 / mineral_mass;
        F_mineral_2_z   = F_mineral_22 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_2.*K_temperature_hi.*dz_hi) ...
                        );              
        mineral_2_conc  = (F_mineral_22./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_2.*K_temperature_hi.*dz_hi) ...
                        );                                                          % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral2 = k_diss_mineral_2.* K_temperature_hi.* mineral_2_conc;     % mineral dissolution rate [mol/m3*y]
        R_diss_int_2     = cumsum(R_diss_mineral2.*dz_hi);                          % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_2     = mineral_2_conc(:)/mineral_2_conc(1);                     % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        % 
        % -----------------------------------------------------------------
        % --- size class 3 ------------------------------------------------
        % ----------------------------------------------------------------- 
        k_sed_min_3 = SSA_3.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                    %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_3 = k_sed_min_3.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));   
        F_mineral_33  =  Flux_mineral_3 / mineral_mass;
        F_mineral_3_z   = F_mineral_33 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_3.*K_temperature_hi.*dz_hi) ...
                        );    
        mineral_3_conc  = (F_mineral_33./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_3.*K_temperature_hi.*dz_hi) ...
                        );                                                           % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral3 = k_diss_mineral_3.* K_temperature_hi.* mineral_3_conc;      % mineral dissolution rate [mol/m3*y]
        R_diss_int_3     = cumsum(R_diss_mineral3.*dz_hi);                           % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_3     = mineral_3_conc(:)/mineral_3_conc(1);                      % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 4 ------------------------------------------------
        % ----------------------------------------------------------------- 
        k_sed_min_4 = SSA_4.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                      %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_4 = k_sed_min_4.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));   
        F_mineral_44  =  Flux_mineral_4 / mineral_mass;
        F_mineral_4_z   = F_mineral_44 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_4.*K_temperature_hi.*dz_hi) ...
                        );                 
        mineral_4_conc  = (F_mineral_44./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_4.*K_temperature_hi.*dz_hi) ...
                        );                                                            % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral4 = k_diss_mineral_4.* K_temperature_hi.* mineral_4_conc;       % mineral dissolution rate [mol/m3*y]
        R_diss_int_4     = cumsum(R_diss_mineral4.*dz_hi);                            % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_4     = mineral_4_conc(:)/mineral_4_conc(1);                       % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 5 ------------------------------------------------
        % -----------------------------------------------------------------
        k_sed_min_5 = SSA_5.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                      %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_5 = k_sed_min_5.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));  
        F_mineral_55  =  Flux_mineral_5 / mineral_mass;
        F_mineral_5_z   = F_mineral_55 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_5.*K_temperature_hi.*dz_hi) ...
                        );               
        mineral_5_conc  = (F_mineral_55./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_5.*K_temperature_hi.*dz_hi) ...
                        );                                                            % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral5 = k_diss_mineral_5.* K_temperature_hi.* mineral_5_conc;       % mineral dissolution rate [mol/m3*y]
        R_diss_int_5     = cumsum(R_diss_mineral5.*dz_hi);                            % depth-integrated dissolution rate [mol/m2*yr]
        % --------------------------- Burial Efficiency -------------------
        BE_mineral_5     = mineral_5_conc(:)/mineral_5_conc(1);                       % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 6 ------------------------------------------------
        % -----------------------------------------------------------------   
        k_sed_min_6 = SSA_6.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                      %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_6 = k_sed_min_6.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));    
        F_mineral_66  =  Flux_mineral_6 / mineral_mass;
        F_mineral_6_z   = F_mineral_66 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_6.*K_temperature_hi.*dz_hi) ...
                        );        
        mineral_6_conc  = (F_mineral_66./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_6.*K_temperature_hi.*dz_hi) ...
                        );                                                            % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral6 = k_diss_mineral_6.* K_temperature_hi.* mineral_6_conc;       % mineral dissolution rate [mol/m3*y]
        R_diss_int_6     = cumsum(R_diss_mineral6.*dz_hi);                            % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_6     = mineral_6_conc(:)/mineral_6_conc(1);                       % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 7 ------------------------------------------------
        % -----------------------------------------------------------------
        k_sed_min_7 = SSA_7.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                      %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_7 = k_sed_min_7.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));  
        F_mineral_77  =  Flux_mineral_7 / mineral_mass;
        F_mineral_7_z   = F_mineral_77 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_7.*K_temperature_hi.*dz_hi) ...
                        );                        
        mineral_7_conc  = (F_mineral_77./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_7.*K_temperature_hi.*dz_hi) ...
                        );                                                            % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral7 = k_diss_mineral_7.* K_temperature_hi.* mineral_7_conc;       % mineral dissolution rate [mol/m3*y]
        R_diss_int_7     = cumsum(R_diss_mineral7.*dz_hi);                            % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_7     = mineral_7_conc(:)/mineral_7_conc(1);                       % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 8 ------------------------------------------------
        % ----------------------------------------------------------------- 
        k_sed_min_8 = SSA_8.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                      %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_8 = k_sed_min_8.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));    
        F_mineral_88  =  Flux_mineral_8 / mineral_mass;
        F_mineral_8_z   = F_mineral_88 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_8.*K_temperature_hi.*dz_hi) ...
                        );                  
        mineral_8_conc  = (F_mineral_88./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_8.*K_temperature_hi.*dz_hi) ...
                        );                                                            % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral8 = k_diss_mineral_8.* K_temperature_hi.* mineral_8_conc;       % mineral dissolution rate [mol/m3*y]
        R_diss_int_8     = cumsum(R_diss_mineral8.*dz_hi);                            % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_8     = mineral_8_conc(:)/mineral_8_conc(1);                       % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 9 ------------------------------------------------
        % -----------------------------------------------------------------         
        k_sed_min_9 = SSA_9.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                      %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_9 = k_sed_min_9.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));    
        F_mineral_99  =  Flux_mineral_9 / mineral_mass;
        F_mineral_9_z   = F_mineral_99 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_9.*K_temperature_hi.*dz_hi) ...
                        );        
        mineral_9_conc  = (F_mineral_99./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_9.*K_temperature_hi.*dz_hi) ...
                        );                                                            % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral9 = k_diss_mineral_9.* K_temperature_hi.* mineral_9_conc;       % mineral dissolution rate [mol/m3*y]
        R_diss_int_9     = cumsum(R_diss_mineral9.*dz_hi);                            % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_9     = mineral_9_conc(:)/mineral_9_conc(1);                       % mineral burial efficiency at depth        
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 10 -----------------------------------------------
        % -----------------------------------------------------------------         
        k_sed_min_10 = SSA_10.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                    %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_10 = k_sed_min_10.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));   
        F_mineral_101  =  Flux_mineral_10 / mineral_mass;
        F_mineral_10_z   = F_mineral_101 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_10.*K_temperature_hi.*dz_hi) ...
                        );                 
        mineral_10_conc  = (F_mineral_101./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_10.*K_temperature_hi.*dz_hi) ...
                        );                                                            % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral10 = k_diss_mineral_10.* K_temperature_hi.* mineral_10_conc;    % mineral dissolution rate [mol/m3*y]
        R_diss_int_10     = cumsum(R_diss_mineral10.*dz_hi);                          % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_10     = mineral_10_conc(:)/mineral_10_conc(1);                    % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 11 -----------------------------------------------
        % ----------------------------------------------------------------- 
        k_sed_min_11 = SSA_11.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                    %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_11 = k_sed_min_11.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));    
        F_mineral_111  =  Flux_mineral_11 / mineral_mass;
        F_mineral_11_z   = F_mineral_111 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_11.*K_temperature_hi.*dz_hi) ...
                        );                       
        mineral_11_conc  = (F_mineral_111./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_11.*K_temperature_hi.*dz_hi) ...
                        );                                                            % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral11 = k_diss_mineral_11.* K_temperature_hi.* mineral_11_conc;    % mineral dissolution rate [mol/m3*y]
        R_diss_int_11     = cumsum(R_diss_mineral11.*dz_hi);                          % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_11     = mineral_11_conc(:)/mineral_11_conc(1);                    % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 12 -----------------------------------------------
        % -----------------------------------------------------------------         
        k_sed_min_12 = SSA_12.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                    %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_12 = k_sed_min_12.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));   
        F_mineral_121  =  Flux_mineral_12 / mineral_mass;
        F_mineral_12_z   = F_mineral_121 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_12.*K_temperature_hi.*dz_hi) ...
                        );              
        mineral_12_conc  = (F_mineral_121./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_12.*K_temperature_hi.*dz_hi) ...
                        );                                                             % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral12 = k_diss_mineral_12.* K_temperature_hi.* mineral_12_conc;     % mineral dissolution rate [mol/m3*y]
        R_diss_int_12     = cumsum(R_diss_mineral12.*dz_hi);                           % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_12     = mineral_12_conc(:)/mineral_12_conc(1);                     % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 13 -----------------------------------------------
        % -----------------------------------------------------------------         
        k_sed_min_13 = SSA_13.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                     %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_13 = k_sed_min_13.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));    
        F_mineral_131  =  Flux_mineral_13 / mineral_mass;
        F_mineral_13_z   = F_mineral_131 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_13.*K_temperature_hi.*dz_hi) ...
                        );                        
        mineral_13_conc  = (F_mineral_131./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_13.*K_temperature_hi.*dz_hi) ...
                        );                                                             % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral13 = k_diss_mineral_13.* K_temperature_hi.* mineral_13_conc;     % mineral dissolution rate [mol/m3*y]
        R_diss_int_13     = cumsum(R_diss_mineral13.*dz_hi);                           % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_13     = mineral_13_conc(:)/mineral_13_conc(1);                     % mineral burial efficiency at depth
        % -----------------------------------------------------------------
        %
        % -----------------------------------------------------------------
        % --- size class 14 -----------------------------------------------
        % -----------------------------------------------------------------         
        k_sed_min_14 = SSA_14.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                     %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral_14 = k_sed_min_14.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));   

        F_mineral_141  =  Flux_mineral_14 / mineral_mass;
        F_mineral_14_z   = F_mineral_141 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_14.*K_temperature_hi.*dz_hi) ...
                        );    
                    
        mineral_14_conc  = (F_mineral_141./(mean(velocity_hi).*365)) ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral_14.*K_temperature_hi.*dz_hi) ...
                        );                                                              % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------
        R_diss_mineral14 = k_diss_mineral_14.* K_temperature_hi.* mineral_14_conc;      % mineral dissolution rate [mol/m3*y]
        R_diss_int_14     = cumsum(R_diss_mineral14.*dz_hi);                            % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------
        BE_mineral_14     = mineral_14_conc(:)/mineral_14_conc(1);                      % mineral burial efficiency at depth
        %
        % -----------------------------------------------------------------
        % --- across all size classes -------------------------------------
        % -----------------------------------------------------------------
        
        mineral_conc   = mineral_1_conc + mineral_2_conc + mineral_3_conc + mineral_4_conc + ...
                          mineral_5_conc + mineral_6_conc + mineral_7_conc + mineral_8_conc + ...
                          mineral_9_conc + mineral_10_conc + mineral_11_conc + mineral_12_conc + ...
                          mineral_13_conc + mineral_14_conc;
        % -----------------------------------------------------------------
        R_diss_mineral = R_diss_mineral1 + R_diss_mineral2 + R_diss_mineral3 + R_diss_mineral4 + ...
                         R_diss_mineral5 + R_diss_mineral6 + R_diss_mineral7 + R_diss_mineral8 + ...
                         R_diss_mineral9 + R_diss_mineral10 + R_diss_mineral11 + R_diss_mineral12 + ...
                         R_diss_mineral13 + R_diss_mineral14;
        % -----------------------------------------------------------------  
        R_diss_int     = cumsum(R_diss_mineral.*dz_hi);                                 % total depth-integrated dissolution rate [mol/m3*y]
        % -----------------------------------------------------------------     
        BE_mineral     = mineral_conc(:)./mineral_conc(1);                              % mineral burial efficiency at depth
        % -----------------------------------------------------------------    
    % ---------------------------------------------------------------------
    else
    % ---------------------------------------------------------------------
        k_sed_min        = SSA.*(60.*60.*24.*365.*(10.^k_min_pH0)).*mineral_mass;                     %1/year ''rate of mineral dissolution''
        alpha_k_mineral  = 10.^(-0.5*log10(org_age_interp) - 0.01);
        k_diss_mineral   = k_sed_min.* exp((-1./(mean(velocity_hi).*365)).*cumsum(alpha_k_mineral.*dz_hi));    
        % -----------------------------------------------------------------
        F_mineral_01  =  Flux_mineral / mineral_mass;
        F_mineral_01_z   = F_mineral_01 ...                                               
                    .*exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral.*K_temperature_hi.*dz_hi) ...
                        );    
                                             
        mineral_conc  = (F_mineral_01./(mean(velocity_hi).*365)) ...                                               
                    *exp((-1./(mean(velocity_hi).*365)).*cumsum(k_diss_mineral.*K_temperature_hi.*dz_hi) ...
                        );                                                              % mineral concentration at depth [mol/m3]
        % -----------------------------------------------------------------    
        R_diss_mineral = k_diss_mineral.* K_temperature_hi.* mineral_conc;              % mineral dissolution rate [mol/m3*y]
        R_diss_int     = cumsum(R_diss_mineral.*dz_hi);                                 % depth-integrated dissolution rate [mol/m2*yr]
        % -----------------------------------------------------------------    
        BE_mineral     = mineral_conc(:)./mineral_conc(1);                              % mineral burial efficiency at depth
        % -----------------------------------------------------------------    
    % ---------------------------------------------------------------------
    end
    % ---------------------------------------------------------------------
    %
    % ---------------------------------------------------------------------
    % --- evaluate convergence criterion for carbonate system -------------
    % ---------------------------------------------------------------------
    K_converge1 = ( mean(CO3_store(iteration_step1,:)) ...
                  -mean(CO3_store(iteration_step1-1,:)) )./...
                  (mean(CO3_store(iteration_step1-1,:)));
    % ---------------------------------------------------------------------
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------
%
% #########################################################################
% ### END OF PARTICLE FACTORY ITERATION ###################################
% #########################################################################
% -------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--- evaluate mineral mass balance ----------------------------------------
%--------------------------------------------------------------------------
rain         = F_mineral;                            % input flux [g/m2/y]
diss         = R_diss_int(end)*mineral_mass;         % dissolution flux [g/m2/y]
bury         = BE_mineral(end)*rain;                 % burial flux [g/m2/y]
% -------------------------------------------------------------------------
mass_balance = ((rain - (diss+bury))/rain)*100;      % mineral mass balance [o/o]
% -------------------------------------------------------------------------
if abs(mass_balance) < 0.1
 disp('*** MASS BALANCE ACHIEVED. ***');
 fprintf(1,'*** mass imbalance = %g percent ***\n',abs(mass_balance));
else
 disp('*** DOUBLE CHECK MASS BALANCE! ***');
 fprintf(1,'*** mass imbalance = %g percent ***\n',abs(mass_balance));
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- MAKE SOME DIAGNOSTIC PLOTS ------------------------------------------
% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
figure(1);
% -------------------------------------------------------------------------
subplot(1,3,1);
plot(T_seawater,z,'k','LineWidth',1.5);
axis ij;
xlim([0 40]);
xlabel('temperature [oC]');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
subplot(1,3,2);
plot(rho_seawater,z,'k','LineWidth',1.5);
axis ij;
xlim([1.025 1.026]);
xlabel('density [g cm^{-3}]');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
subplot(1,3,3);
plot(nho,z,'k','LineWidth',1.5);
axis ij;
%xlim([0.00 0.02]);
xlabel('kinematic viscosity [m^2 d^{-1}]');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
figure(2);
% -------------------------------------------------------------------------
subplot(1,2,1);
plot(velocity_total,z,'k','LineWidth',1.5);
axis ij;
xlabel('mean settling velocity [m day^{-1}]');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
subplot(1,2,2);
plot(org_age_d,z,'k','LineWidth',1.5);
axis ij;
xlabel('C_{org} age [day]');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
Z0_data_Store = [z' velocity_total' org_age_d'];
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
figure(3);
% -------------------------------------------------------------------------
subplot(2,2,1);
plot(Rate_org_int,z_hi,'k','LineWidth',1.5);
axis ij;
xlabel('C_{org} remineralization [mol m^{-2} y^{-1}]');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
subplot(2,2,2);
plot(BE_org,z_hi,'k','LineWidth',1.5);
axis ij;
xlabel('C_{org} burial efficiency');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
subplot(2,2,3);
plot(R_diss_int,z_hi,'k','LineWidth',1.5);
axis ij;
xlabel('integrated dissolution [mol m^{-2} y^{-1}]');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
subplot(2,2,4);
plot(BE_mineral,z_hi,'k','LineWidth',1.5);
axis ij;
xlabel('silicate burial efficiency');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
Z1_data_Store = [z_hi' Rate_org_int' BE_org R_diss_int' BE_mineral];
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
figure(4);
% -------------------------------------------------------------------------
n_plot1 = 6;
% -------------------------------------------------------------------------
subplot(1,n_plot1,1)
plot(Oxygen, z,'k','LineWidth',1.5); axis ij
xlabel('[O_2] (\muM)');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
if strcmp(Calibration_site,'GofM A')
% -------------------------------------------------------------------------
subplot(1,n_plot1,2)
plot(C_DIC, z,'k',1000.*DIC_1_Validation(:,1),DIC_1_Validation(:,2),'ob','LineWidth',1.5); axis ij
xlabel('[DIC] (\muM)');
ylabel('depth [m]');
ylim([0 500]);
box on; grid on;
% -------------------------------------------------------------------------
subplot(1,n_plot1,3)
plot(ALK, z,'k',1000.*ALK_1_Validation(:,1),ALK_1_Validation(:,2),'ob','LineWidth',1.5); axis ij 
xlabel('[TA] (\muM)');
ylabel('depth [m]');
ylim([0 500]);
box on; grid on;
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
if strcmp(Calibration_site,'GofM B')
% -------------------------------------------------------------------------
subplot(1,n_plot1,2)
plot(C_DIC, z,'k',1000.*DIC_2_Validation(:,1),DIC_2_Validation(:,2),'ob','LineWidth',1.5); axis ij
xlabel('[DIC] (\muM)');
ylabel('depth [m]');
ylim([0 500]);
box on; grid on;
% -------------------------------------------------------------------------
subplot(1,n_plot1,3)
plot(ALK, z,'k',1000.*ALK_2_Validation(:,1),ALK_2_Validation(:,2),'ob','LineWidth',1.5); axis ij 
xlabel('[TA] (\muM)');
ylabel('depth [m]');
ylim([0 500]);
box on; grid on;
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
if strcmp(Calibration_site,'Pacific')
% -------------------------------------------------------------------------
subplot(1,n_plot1,2)
plot(C_DIC, z,'k','LineWidth',1.5); axis ij
xlabel('[DIC] (\muM)');
ylabel('depth [m]');
ylim([0 500]);
box on; grid on;
% -------------------------------------------------------------------------
subplot(1,n_plot1,3)
plot(ALK, z,'k',Cali_ALK(:,1),Cali_ALK(:,2),'ob','LineWidth',1.5); axis ij 
xlabel('[TA] (\muM)');
ylabel('depth [m]');
ylim([0 500]);
box on; grid on;
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
if strcmp(Calibration_site,'Pacific')
 % -------------------------------------------------------------------------   
subplot(1,n_plot1,4)   
plot(pH, z,'k',Cali_pH(:,1),Cali_pH(:,2),'ob','LineWidth',1.5); axis ij
xlabel('[pH] (\muM)');
ylabel('depth [m]');
ylim([0 500]);
box on; grid on;
% -------------------------------------------------------------------------
else
   % ------------------------------------------------------------------------- 
subplot(1,n_plot1,4)
plot(pH, z,'k','LineWidth',1.5); axis ij
xlabel('[pH] (\muM)');
ylabel('depth [m]');
box on; grid on;
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
subplot(1,n_plot1,5)
plot(CO3_1, z,'k','LineWidth',1.5); axis ij
xlabel('[CO_3] (\muM)');
ylabel('depth [m]'); 
box on; grid on;
% -------------------------------------------------------------------------
subplot(1,n_plot1,6)
plot(sigma_carb, z,'k','LineWidth',1.5); axis ij
xlabel('\Omega-1');
ylabel('depth [m]'); 
box on; grid on;
% -------------------------------------------------------------------------
%
% ------------------------------------------------------------------------- 
Z2_data_Store = [z' C_DIC' ALK' pH' CO3_1' sigma_carb'];
% -------------------------------------------------------------------------
AAA_time = toc;
% -------------------------------------------------------------------------
%