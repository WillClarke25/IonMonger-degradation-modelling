function params = paper_params(rec_fac,FREQ,EXP,flag)
% This function returns a parameter structure that contains all the inputs
% required for the simulation. rec_fac is a function handle determining the
% time-dependent recombination factor associated with degradation. This
% parameter structure computes a single frequency measurement.

% Defaults to recombination parameter set A
if nargin<4, flag='A'; end

%% Program settings

% Options for the user
workfolder  = './Data/'; % the folder to which data will be saved,
% note that the string must end with a forward slash
OutputFcn   = 'PrintTime'; % ode15s optional message function, choose
% either 'PrintVolt' or 'PrintTime' (which can be found in Code/Common/)
Verbose     = false; % set this to false to suppress message output,
% note that this option overwrites the previous two if Verbose=false
UseSplits   = true; % set this to false to make a single call to ode15s

% Resolution and error tolerances
N    = 400; % Number of subintervals, with N+1 being the number of grid points
rtol = 1e-6; % Relative temporal tolerance for ode15s solver
atol = 1e-6; % Absolute temporal tolerance for ode15s solver
phidisp = 100; % displacement factor for electric potential (VT) (optional)

%% Parameter input

% Physical constants
eps0  = 8.854187817e-12;  % permittivity of free space (Fm-1)
q     = 1.6021766209e-19; % charge on a proton (C)
Fph   = 1.4e21;           % incident photon flux (m-2s-1)
kB    = 8.61733035e-5;    % Boltzmann constant (eVK-1)

% Perovskite parameters
T     = 298;       % temperature (K)
b     = 400e-9;    % perovskite layer width (m) (normally between 150-600nm)
epsp  = 24.1*eps0; % permittivity of perovskite (Fm-1)
alpha = 1.3e7;     % perovskite absorption coefficient (m-1)
Ec    = -3.7;      % conduction band minimum (eV)
Ev    = -5.4;      % valence band maximum (eV)
Dn    = 1.7e-4;    % perovskite electron diffusion coefficient (m2s-1)
Dp    = 1.7e-4;    % perovskite hole diffusion coefficient (m2s-1)
gc    = 8.1e24;    % conduction band density of states (m-3)
gv    = 5.8e24;    % valence band density of states (m-3)

% Ion parameters
N0    = 5e23;        % typical density of ion vacancies (m-3)
D     = @(Dinf, EA) Dinf*exp(-EA/(kB*T)); % diffusivity relation
DI    = 3.5e-14; %D(DIinf, EAI); % diffusion coefficient for iodide ions (m2s-1)

% Direction of light
inverted = false; % choose false for a standard architecture cell (light
% entering through the ETL), true for an inverted architecture cell
% (light entering through the HTL)

% ETL parameters
dE    = 1e24;    % effective doping density of ETL (m-3)
gcE   = 5e25;    % effective conduction band DoS in ETL (m-3)
EcE   = -4.0;    % conduction band reference energy in ETL (eV)
bE    = 100e-9;  % width of ETL (m)
epsE  = 10*eps0; % permittivity of ETL (Fm-1)
DE    = 1e-5;    % electron diffusion coefficient in ETL (m2s-1)

% HTL parameters
dH    = 1e24;    % effective doping density of HTL (m-3) 
gvH   = 5e25;    % effective valence band DoS in HTL (m-3)
EvH   = -5.1;    % valence band reference energy in HTL (eV)
bH    = 200e-9;  % width of HTL (m)
epsH  = 3*eps0;  % permittivity of HTL (Fm-1)
DH    = 1e-6;    % hole diffusion coefficient in HTL (m2s-1)

if strcmp(flag,'A') % Recombination parameter set A
    % Bulk recombination
    tn    = 3e-9;    % electron pseudo-lifetime for SRH (s)
    tp    = 3e-7;    % hole pseudo-lifetime for SRH (s)
    beta  = 1e-15;       % bimolecular recombination rate (m3s-1)
    Augn  = 0;       % electron-dominated Auger recombination rate (m6s-1)
    Augp  = 0;       % hole-dominated Auger recombination rate (m6s-1)
    
    % Interface recombination (max. velocity ~ 1e5)
    betaE = 0;       % effective ETL/perovskite bimolecular rate (m3s-1)
    betaH = 0;       % effective perovskite/HTL bimolecular rate (m3s-1)
    vnE   = 1e5;     % effective electron recombination velocity for SRH (ms-1)
    vpE   = 10;      % hole recombination velocity for SRH (ms-1)
    vnH   = 0.1;     % electron recombination velocity for SRH (ms-1)
    vpH   = 1e5;     % effective hole recombination velocity for SRH (ms-1)
elseif strcmp(flag,'B') % Recombination parameter set A
    % Bulk recombination
    tn    = 1e-10;    % electron pseudo-lifetime for SRH (s)
    tp    = 1e-7;    % hole pseudo-lifetime for SRH (s)
    beta  = 1e-16;       % bimolecular recombination rate (m3s-1)
    Augn  = 0;       % electron-dominated Auger recombination rate (m6s-1)
    Augp  = 0;       % hole-dominated Auger recombination rate (m6s-1)
    
    % Interface recombination (max. velocity ~ 1e5)
    betaE = 0;       % effective ETL/perovskite bimolecular rate (m3s-1)
    betaH = 0;       % effective perovskite/HTL bimolecular rate (m3s-1)
    vnE   = 0e5;     % effective electron recombination velocity for SRH (ms-1)
    vpE   = 0;      % hole recombination velocity for SRH (ms-1)
    vnH   = 0;     % electron recombination velocity for SRH (ms-1)
    vpH   = 0e5;     % effective hole recombination velocity for SRH (ms-1)
end

%% Option to set initial distributions from a saved solution

% Name of a saved solution structure. Initial distributions will be taken
% from the final distributions of the saved solution. In this case,
% `applied_voltage` does not need an initial voltage.

% input_filename = 'Data/simulation.mat';

%% Non-dimensionalise model parameters and save all inputs

% Compile all parameters into a convenient structure
vars = setdiff(who,{'params','vars'});
for i=1:length(vars), params.(vars{i}) = eval(vars{i}); end

% Non-dimensionalise the user-defined input parameters
params = nondimensionalise(params);

% Unpack variables and functions needed in the rest of this function
[tstar2t, psi2Vap, Upsilon, Vbi] = struct2array(params, ...
    {'tstar2t','psi2Vap','Upsilon','Vbi'});

%% Simulation protocol
% In order to make use of construct_protocol.m, instructions must be given
% in a specific order and format. Please see the GUIDE.md. Otherwise one
% can specify their own dimensionless functions of time (light and psi),
% dimensionless vectors (time and spEDRFTG6YH7UJ8IKOLP[;'~
% lits) and option whether to findVoc.

% Light protocol
light_intensity = {1};

% Voltage protocol
nwaves = 8; % number of waves simulated during the integration time
min_cycles = 2; % minimum number of cycles that must ALWAYS be measured
Tint = 5; % integration time (s)
Vp = 10e-3; % AC voltage amplitude (V)
rest_time = 1; % rest time between frequencies (s)

if strcmp(flag,'A') ; VDC = 0.9;
elseif strcmp(flag,'B') ; VDC = 1.05;
end

if strcmp(EXP,'IS')
    % sequential impedance spectrum measurement
    applied_voltage = {VDC};

    if nwaves/FREQ < Tint
        num_cycles = nwaves;
        dummy_time = Tint-nwaves/FREQ;
    elseif min_cycles/FREQ < Tint
        num_cycles = min_cycles;
        dummy_time = Tint-min_cycles/FREQ;
    else
        num_cycles = min_cycles;
        dummy_time = 0;
    end

    % Linear DC section
    applied_voltage{end+1} = 'linear';
    applied_voltage{end+1} = rest_time+dummy_time;
    applied_voltage{end+1} = VDC;

    % sinusoidal sections
    for j = 1:num_cycles
        applied_voltage{end+1} = 'sin';
        applied_voltage{end+1} = 1/FREQ;
        applied_voltage{end+1} = VDC+Vp;
    end
elseif strcmp(EXP,'jv')
    % high-resolution JV sweep
    scan_rate = FREQ; % use FREQ variable as scan rate (V/s)
    t_step = 0.1/scan_rate;
    applied_voltage = {Vbi,'tanh',20,1.2,...
        'linear',t_step,1.1,...
        'linear',t_step,1.0,...
        'linear',t_step,0.9,...
        'linear',t_step,0.8,...
        'linear',t_step,0.7,...
        'linear',t_step,0.6,...
        'linear',t_step,0.5,...
        'linear',t_step,0.4,...
        'linear',t_step,0.3,...
        'linear',t_step,0.2,...
        'linear',t_step,0.1,...
        'linear',t_step,0,...
        'linear',10*t_step,1.2}; 
elseif strcmp(EXP,'parallel')
    % impedance spectroscopy using IonMonger's default parallel approach
    % WARNING: this protocol can only be used when there is no degradation
    % and FREQ becomes the number of frequencies to sample
    applied_voltage = {'impedance', ...
        1e-1, ...     % minimum impedance frequency (Hz)
        1e6, ...      % maximum impedance frequency (Hz)
        VDC, ...      % DC voltage (V)
        Vp, ...    % AC voltage amplitude (V)
        FREQ, ...       % number of frequencies to sample
        5};           % number of sine waves
elseif strcmp(EXP,'fixed_voltage')
    applied_voltage = {VDC,'linear',FREQ,VDC};
end
reduced_output = true; % set to true to reduce the amount of data retained ...
% in impedance simulations

% Choose whether the time points are spaced linearly or logarithmically
time_spacing = 'lin'; % set equal to either 'lin' (default) or 'log'

%% Create the simulation protocol and plot (if Verbose)

% Create the protocol and time points automatically, psi=(Vbi-Vap)/(2*kB*T)
[light, psi, time, splits, findVoc] = ...
    construct_protocol(params,light_intensity,applied_voltage,time_spacing);

if strcmp(EXP,'IS')
    % remove unnecessary split points
    split_indices = [1];
    split_indices(end+1) = split_indices(end)+1;
    split_indices(end+1) = split_indices(end)+num_cycles;
    splits = splits(split_indices);
end

% Apply the options defined above
if inverted, inv = -1; else, inv = 1; end
if ~UseSplits, splits = time([1,end]); end

% Define the charge carrier generation function G(x,t)
G = @(x,t) light(t).*Upsilon./(1-exp(-Upsilon)).*exp(-Upsilon*(inv*x+(1-inv)/2));

% Plot the light regime
if Verbose
    if ishandle(98), clf(98); end; figure(98);
    plot(tstar2t(time),light(time),'Color',[0.93 0.69, 0.13]);
    xlabel('Time (s)'); ylabel('Light intensity (Sun equiv.)');
    title('light(t)');
    drawnow;
end

% Plot the voltage regime
if Verbose
    if ishandle(99), clf(99); end; figure(99);
    if isnan(psi(time(end)))
        title('Open-circuit');
    else
        plot(tstar2t(time),psi2Vap(psi(time)));
        xlabel('Time (s)'); ylabel('Applied Voltage (V)');
        title('V(t)');
        if findVoc
            hold on; plot(0,Vbi,'o','MarkerSize',8);
            title('V(t) except the voltage starts from Voc, not Vbi as shown here');
        end
        if strcmp(applied_voltage{1},'impedance')
            title({'example impedance protocol','at a frequency of 1Hz'});
            hold on;
            plot(tstar2t(time(end-200:end)),psi2Vap(psi(time(end-200:end))),'r');
            legend({'','phase analysis region'},'Location','best');
            hold off;
        end
        ylim([min([0,psi2Vap(psi(time))]), ceil(2*max(psi2Vap(psi(time))))/2]);
    end
    drawnow;
end


%% Compile more parameters into the params structure
vars = setdiff(setdiff(who,fieldnames(params)),{'params','vars','i'});
for i=1:length(vars), params.(vars{i}) = eval(vars{i}); end

% Make the folder in which to save the output data (specified above)
if exist(workfolder,'dir')~=7, mkdir(workfolder); end


end
