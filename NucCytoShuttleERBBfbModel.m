function [OutputTable,CellVolumeTable] = NucCytoShuttleERBBfbModel(Spike,varargin)
%NucCytoShuttleModel    Simulates a spike of reporter cargo in a cell type
%                       with known initial conditions and nuclear-cytoplasmic 
%                       volumes.
%
%INPUTS
%    Spike: Concentration of induced or microinjected cargo in µM.  No 
%        default value.
%    CellType: Two-element vector of [Nuclear Cytoplasmic] volumes in pL.
%        Default = [0.45; 1.55] (HeLa estimates).
%    CellConcs: Ten-element vector of total cellular concentrations of [Ran
%        RanBP1 RanGap RCC1 NTF2 Impα Impβ Cas CRM1 RanBP3] in µM.
%        Default = [5; 2; 0.5; 0.25; 0.6; 1; 3; 3; 0.3; 0.05] (HeLa estimates).
%    LumpedReceptors = Four-element vector of approximate total cellular
%        concentrations of [α-β cargo; β cargo; other transport receptors;
%        NES cargo] in µM. Default = [10; 1; 3.6; 1] (HeLa estimates).
%    NUPs = Nuclear pores in copies per cell. Empty vector indicates
%        modeling of NUP-unlimited transport, scalar indicates modeling of
%        cargo-universal NUPs, two-element vector indicates modeling of
%        [cargo-universal NUPs cargo-excluded NUPs].  Default = 3000 (HeLa
%        estimates).
%    NLSparams = NLS kinetic parameters for sweeping through with ErbB
%        feedback.  Four element vector indicates NLS paramters for
%        [kon1_ABCBP80; koff1_ABCBP80; kon2_ABCBP80; koff2_ABCBP80].  
%        Default = [0.15; 0.01; 0.0064; 0.00025];
%
%OUTPUTS
%    OutputTable: Table of seven tables corresponding to the mass-per-cell
%        spike-in trajectories for NLS, IBB, CBP80, NNES, INES, and CNES 
%        cargo, plus the initial trajectory to steady state (SS).
%    CellVolumeTable: Table of three entries for the calculated cytoplasmic
%        volume (VolCyto) and nuclear volume (VolNuc) after accounting for
%        the NPC compartment (VolNPC).

%% Set defaults and parse inputs

%Set input parameter defaults
defaultCellType = [0.45; 1.55]; %pL from Jarnac code of Riddick & Macara, J Cell Biol (2005)
defaultCellConcs = [5; 2; 0.5; 0.25; 0.6; 1; 3; 3; 0.3; 0.05]; %µM from Table I of Riddick & Macara, J Cell Biol (2005) and Quantitative Western Blotting text of Riddick & Macara, Mol Syst Biol (2007)
defaultLumpedReceptors = [10; 1; 3.6; 1]; %µM from Table I of Riddick & Macara, J Cell Biol (2005) and cross-referencing of la Cour et al., Nucl Acids Res (2003) with Liu et al., Nat Biotechnol (2019)
defaultNUPs = 3000; %copies per cell Maul et al., J Cell Biol (1972)
defaultNLSparams = [0.15; 0.01; 0.0064; 0.00025];

%Set criteria for inputs
validSpike = @(x) isnumeric(x) & (x >= 0); %Spike should be ≥ 0 µM
validCellType = @(x) ((isrow(x) & size(x,2) == 2) | (iscolumn(x) ...
    & size(x,1) == 2)) & all(x > 0); %CellType should be a two-element row or column vector of volumes > 0 pL
validCellConcs = @(x) ((isrow(x) & size(x,2) == 10) | (iscolumn(x) ...
    & size(x,1) == 10)) & all(x >= 0); %CellConcs should be a ten-element row or column vector of volumes ≥ 0 µM
validLumpedReceptors = @(x) ((isrow(x) & size(x,2) == 4) | (iscolumn(x) ...
    & size(x,1) == 4)) & all(x >= 0); %LumpedParameters should be a four-element row or column vector of volumes ≥ 0 µM
validNUPs = @(x) all(isempty(x)) | all((isscalar(x) & x >= 0)) | (((isrow(x) & ...
    size(x,2) == 2) | (iscolumn(x) & size(x,1) == 2)) & all(x >= 0)); %NUPs should be empty, a scalar ≥ 0 copies per cell, or a two-element row or column vector of ≥ copies per cell
validNLSparams = @(x) ((isrow(x) & size(x,2) == 4) | (iscolumn(x) ...
    & size(x,1) == 4)) & all(x >= 0); %NLSparams should be a four-element row or column vector of paramters ≥ 0

%Set parser, parser options, and valid input criteria
NucCytoShuttleModelparser = inputParser;
addRequired(NucCytoShuttleModelparser,'Spike',validSpike);
addOptional(NucCytoShuttleModelparser,'CellType',defaultCellType,validCellType);
addOptional(NucCytoShuttleModelparser,'CellConcs',defaultCellConcs,validCellConcs);
addOptional(NucCytoShuttleModelparser,'LumpedReceptors',defaultLumpedReceptors,validLumpedReceptors);
addOptional(NucCytoShuttleModelparser,'NUPs',defaultNUPs,validNUPs);
addOptional(NucCytoShuttleModelparser,'NLSparams',defaultNLSparams,validNLSparams);

%Parse and extract inputs
parse(NucCytoShuttleModelparser,Spike,varargin{:});
Spike = NucCytoShuttleModelparser.Results.Spike;
CellType = NucCytoShuttleModelparser.Results.CellType;
CellConcs = NucCytoShuttleModelparser.Results.CellConcs;
LumpedReceptors = NucCytoShuttleModelparser.Results.LumpedReceptors;
NUPs = NucCytoShuttleModelparser.Results.NUPs;
NLSparams = NucCytoShuttleModelparser.Results.NLSparams;

%% Convert inputs to column vectors
if isrow(CellType)
    CellType = CellType';
end
if isrow(CellConcs)
    CellConcs = CellConcs';
end
if isrow(LumpedReceptors)
    LumpedReceptors = LumpedReceptors';
end
if isempty(NUPs)
    NUPs = [0; 0];
end
if isscalar(NUPs)
    NUPs = [NUPs 0];
end
if isrow(NUPs)
    NUPs = NUPs';
end

%% Calculate cell compartments, input masses for rate equations, and read in rate parameters

CellVolumeTable = cellcompart(CellType); %Calculate cell compartments
NUPs = NUPs/6.02e23/sum(CellVolumeTable{1,:})*1e18; %Convert NUP copy number to concentration; convert mol/pL to µmol/L (µM)
IMsTable = readtabletranspose(strcat(pwd,'/InitialConditionsDefault.xlsx')); %Read in default initial conditions 
IMsTable{1,1:16} = [CellConcs' LumpedReceptors' NUPs']*sum(CellVolumeTable{1,:}); %Replace concentrations with masses
IMsTable.Properties.VariableUnits(:) = {'amol'}; %Update units
ParameterTable = readtabletranspose(strcat(pwd,'/ModelParametersDefault.xlsx')); %Read in default model parameters; units described in 'VariableUnits'

%% Update NLS parameters according to NLSparams –> assign to CBP80 because NLS parameters are tied to RanBP3 kinetics
ParameterTable.kon1_ABCBP80 = NLSparams(1);
ParameterTable.koff1_ABCBP80 = NLSparams(2);
ParameterTable.kon2_ABCBP80 = NLSparams(3);
ParameterTable.koff2_ABCBP80 = NLSparams(4);

%% Simulate to steady state

NucCytoODEopts = odeset('RelTol',1e-8,'AbsTol',1e-8,'NonNegative',1:width(IMsTable));
tspan = linspace(0,3000,3001); %Simulate for 3000 sec (sufficient for all species to reach steady state)
[~,Yss] = ode15s(@(t,y)NucCytoShuttleERBBfbODE(t,y,IMsTable.Properties.VariableNames, ...
    ParameterTable{1,:},ParameterTable.Properties.VariableNames, ...
    CellVolumeTable{1,:},~all(NUPs == 0),true),tspan,IMsTable{1,:},NucCytoODEopts);
SS = array2table(Yss,'VariableNames',IMsTable.Properties.VariableNames);
SS.Properties.VariableUnits(:) = {'amol'}; %Add units

%% Simulate cargo spike-in to steady state

%CBP80-like cargo
IMs_CBP80 = Yss(height(Yss),:); %Initial masses from steady-state output
IMs_CBP80(strcmp(IMsTable.Properties.VariableNames,'CBP80y')) = Spike*sum(CellVolumeTable{1,:}); %Input mass of spike-in cargo
[~,Y_CBP80] = ode15s(@(t,y)NucCytoShuttleERBBfbODE(t,y,IMsTable.Properties.VariableNames, ...
    ParameterTable{1,:},ParameterTable.Properties.VariableNames, ...
    CellVolumeTable{1,:},~all(NUPs == 0),true),tspan,IMs_CBP80,NucCytoODEopts);
CBP80 = array2table(Y_CBP80,'VariableNames',IMsTable.Properties.VariableNames);
CBP80.Properties.VariableUnits(:) = {'amol'}; %Add units

OutputTable = table(CBP80,SS);
OutputTable.Properties.VariableUnits(:) = {'amol'}; %Add units

%% Local functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t = cellcompart(v)
%cellcompart    Calculates cytoplasmic, nuclear, and NPC volumes
%t = cellcompart(v)
%t: output table of cytoplasmic, nuclear, and NPC volumes in pL
%v: input volumes of NPC-free nuclear and cytoplasmic volumes in pL

NPCheight = 70e-9; %Height of human NPC in m from Figure 2A of Bui et al., Cell (2013)

vnew(3) = (4*pi/3*((v(1)*1e-12*3/4/pi)^(1/3)+NPCheight/2)^3 ...
    -4*pi/3*((v(1)*1e-12*3/4/pi)^(1/3)-NPCheight/2)^3)*1e12; %Calculate annular shell around nucleus with height of the NPC
vnew(1) = v(1) + v(2) - (4*pi/3*((v(1)*1e-12*3/4/pi)^(1/3)+NPCheight/2)^3)*1e12; %Update cytoplasmic volume
vnew(2) = (4*pi/3*((v(1)*1e-12*3/4/pi)^(1/3)-NPCheight/2)^3)*1e12; %Update nuclear volume

t = array2table(vnew,'VariableNames',{'VolCyto' 'VolNuc' 'VolNPC'});
t.Properties.VariableUnits(:) = {'pL'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t = readtabletranspose(fn)
%Reads in table and transposes rows-columns as variable names.  Expected
%format is the variable name (Column A), number (Column B), units
%(Column C), and notes as needed (Column D)
t_temp = readtable(fn, 'ReadRowNames', 1);
if width(t_temp) == 1
    t = array2table(t_temp.Variables','VariableNames',t_temp.Properties.RowNames');
elseif width(t_temp) == 2
    t = array2table(t_temp{:,1}','VariableNames',t_temp.Properties.RowNames');
    t.Properties.VariableUnits = t_temp{:,2}'; %Extract units
else
    t = array2table(t_temp{:,1}','VariableNames',t_temp.Properties.RowNames');
    t.Properties.VariableUnits = t_temp{:,2}'; %Extract units
    t.Properties.VariableDescriptions = t_temp{:,3}'; %Extract notes and discard anything else
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%


end