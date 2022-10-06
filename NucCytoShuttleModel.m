function [OutputTable,CellVolumeTable] = NucCytoShuttleModel(Spike,varargin)
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

%Set parser, parser options, and valid input criteria
NucCytoShuttleModelparser = inputParser;
addRequired(NucCytoShuttleModelparser,'Spike',validSpike);
addOptional(NucCytoShuttleModelparser,'CellType',defaultCellType,validCellType);
addOptional(NucCytoShuttleModelparser,'CellConcs',defaultCellConcs,validCellConcs);
addOptional(NucCytoShuttleModelparser,'LumpedReceptors',defaultLumpedReceptors,validLumpedReceptors);
addOptional(NucCytoShuttleModelparser,'NUPs',defaultNUPs,validNUPs);

%Parse and extract inputs
parse(NucCytoShuttleModelparser,Spike,varargin{:});
Spike = NucCytoShuttleModelparser.Results.Spike;
CellType = NucCytoShuttleModelparser.Results.CellType;
CellConcs = NucCytoShuttleModelparser.Results.CellConcs;
LumpedReceptors = NucCytoShuttleModelparser.Results.LumpedReceptors;
NUPs = NucCytoShuttleModelparser.Results.NUPs;

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
% IMsTable = readtabletranspose(strcat(pwd,'/InitialConditionsDefault.xlsx')); %Read in default initial conditions 
IMsTable = array2table(zeros(1,207),'VariableNames',{'Rd' 'P' 'RanGap' ...
    'RCC1' 'Ny' 'Ay' 'By' 'Sy' 'My' 'P3y' 'Cy' 'C2y' 'NRy' 'NESy' 'NUP1' ...
    'NUP2' 'Rt' 'RdNn' 'Rdn' 'Nn' 'C2n' 'Sn' 'Cn' 'An' 'Bn' 'ABn' 'ABCn' ...
    'RtBn' 'RtSAn' 'C2Bn' 'ABCFn' 'ACFn' 'CFn' 'C2Fn' 'C2FBn' 'RCC1Rd' ...
    'RCC1R' 'RCC1Rt' 'ABntemp' 'NRn' 'RtNRn' 'ACn' 'RtSn' 'RdNy' 'ABy' ...
    'ABCy' 'RtBy' 'RtSAy' 'C2By' 'RtPBy' 'CFy' 'ABCFy' 'C2Fy' 'C2FBy' ...
    'ABCytemp' 'ABCFytemp' 'C2Bytemp' 'C2FBytemp' 'RtNRy' 'Rty' 'RtPy' ...
    'ABytemp' 'RtNRPy' 'CBP80y' 'NNESy' 'INESy' 'CNESy' 'ABCBP80ytemp' ...
    'ABNNESytemp' 'ABCNESytemp' 'ABCBP80y' 'ABCNESy' 'BINESy' 'BINESytemp' ...
    'ABNNESy' 'XPO1y' 'XPO1Ny' 'XPO1Iy' 'XPO1Cy' 'CBP80n' 'NESn' 'NNESn' ...
    'INESn' 'CNESn' 'Mn' 'P3n' 'MP3n' 'RtMn' 'RtP3n' 'RtMP3n' 'MP3NESn' ...
    'MP3NNESn' 'MP3INESn' 'MP3CNESn' 'XPO1n' 'XPO1Cn' 'XPO1In' 'XPO1Nn' ...
    'ACBP80n' 'ACNESn' 'BINESn' 'ANNESn' 'ABCBP80n' 'ABCNESn' 'ABNNESn' ...
    'Byperi' 'Nyperi' 'Syperi' 'Myperi' 'C2Byperi' 'C2FByperi' 'RdNyperi' ...
    'RtByperi' 'RtSAyperi' 'AByperi' 'ABCyperi' 'ABCFyperi' 'ABCBP80yperi' ...
    'ABNNESyperi' 'ABCNESyperi' 'BINESyperi' 'XPO1Nyperi' 'XPO1Iyperi' ...
    'XPO1Cyperi' 'Bnperi' 'Nnperi' 'Snperi' 'Mnperi' 'C2Bnperi' 'C2FBnperi' ...
    'RdNnperi' 'RtBnperi' 'RtSAnperi' 'ABnperi' 'ABCnperi' 'ABCFnperi' ...
    'ABCBP80nperi' 'ABNNESnperi' 'ABCNESnperi' 'BINESnperi' 'XPO1Nnperi' ...
    'XPO1Inperi' 'XPO1Cnperi' 'NUP1_Bn' 'NUP1_Nn' 'NUP1_Sn' 'NUP1_Mn' ...
    'NUP1_C2Bn' 'NUP1_C2FBn' 'NUP1_RdNn' 'NUP1_RtBn' 'NUP1_RtSAn' ...
    'NUP1_ABn' 'NUP1_ABCn' 'NUP1_ABCFn' 'NUP1_ABCBP80n' 'NUP1_ABNNESn' ...
    'NUP1_ABCNESn' 'NUP1_BINESn' 'NUP1_XPO1Nn' 'NUP1_XPO1In' 'NUP1_XPO1Cn' ...
    'NUP1_By' 'NUP1_Ny' 'NUP1_Sy' 'NUP1_My' 'NUP1_C2By' 'NUP1_C2FBy' ...
    'NUP1_RdNy' 'NUP1_RtBy' 'NUP1_RtSAy' 'NUP1_ABy' 'NUP1_ABCy' 'NUP1_ABCFy' ...
    'NUP1_ABCBP80y' 'NUP1_ABNNESy' 'NUP1_ABCNESy' 'NUP1_BINESy' 'NUP1_XPO1Ny' ...
    'NUP1_XPO1Iy' 'NUP1_XPO1Cy' 'ABP3y' 'ABP3ytemp' 'ABP3yperi' 'NUP1_ABP3y' ...
    'ABP3n' 'ABP3nperi' 'NUP1_ABP3n' 'AP3n' 'NUP2_Bn' 'NUP2_Nn' 'NUP2_Sn' ...
    'NUP2_Mn' 'NUP2_RdNn' 'NUP2_RtBn' 'NUP2_RtSAn' 'NUP2_ABn' 'NUP2_By' ...
    'NUP2_Ny' 'NUP2_Sy' 'NUP2_My' 'NUP2_RdNy' 'NUP2_RtBy' 'NUP2_RtSAy' ...
    'NUP2_ABy' 'NUP2_ABP3y' 'NUP2_ABP3n'}); %Initialize initial conditions table
IMsTable{1,1:16} = [CellConcs' LumpedReceptors' NUPs']*sum(CellVolumeTable{1,:}); %Replace concentrations with masses
IMsTable.Properties.VariableUnits(:) = {'amol'}; %Update units
ParameterTable = readtabletranspose(strcat(pwd,'/ModelParametersDefault.xlsx')); %Read in default model parameters; units described in 'VariableUnits'

%% Simulate to steady state

NucCytoODEopts = odeset('RelTol',1e-8,'AbsTol',1e-8,'NonNegative',1:width(IMsTable));
tspan = linspace(0,3000,3001); %Simulate for 3000 sec (sufficient for all species to reach steady state)
[~,Yss] = ode15s(@(t,y)NucCytoShuttleODE(t,y,IMsTable.Properties.VariableNames, ...
    ParameterTable{1,:},ParameterTable.Properties.VariableNames, ...
    CellVolumeTable{1,:},~all(NUPs == 0)),tspan,IMsTable{1,:},NucCytoODEopts);
SS = array2table(Yss,'VariableNames',IMsTable.Properties.VariableNames);
SS.Properties.VariableUnits(:) = {'amol'}; %Add units

%% Simulate cargo spike-ins to steady state

%NLS-like cargo
IMs_CF = Yss(height(Yss),:); %Initial masses from steady-state output
IMs_CF(strcmp(IMsTable.Properties.VariableNames,'CFy')) = Spike*sum(CellVolumeTable{1,:}); %Input mass of spike-in cargo
[~,Y_CF] = ode15s(@(t,y)NucCytoShuttleODE(t,y,IMsTable.Properties.VariableNames, ...
    ParameterTable{1,:},ParameterTable.Properties.VariableNames, ...
    CellVolumeTable{1,:},~all(NUPs == 0)),tspan,IMs_CF,NucCytoODEopts);
NLS = array2table(Y_CF,'VariableNames',IMsTable.Properties.VariableNames);
NLS.Properties.VariableUnits(:) = {'amol'}; %Add units

%IBB-like cargo
IMs_C2F = Yss(height(Yss),:); %Initial masses from steady-state output
IMs_C2F(strcmp(IMsTable.Properties.VariableNames,'C2Fy')) = Spike*sum(CellVolumeTable{1,:}); %Input mass of spike-in cargo
[~,Y_C2F] = ode15s(@(t,y)NucCytoShuttleODE(t,y,IMsTable.Properties.VariableNames, ...
    ParameterTable{1,:},ParameterTable.Properties.VariableNames, ...
    CellVolumeTable{1,:},~all(NUPs == 0)),tspan,IMs_C2F,NucCytoODEopts);
IBB = array2table(Y_C2F,'VariableNames',IMsTable.Properties.VariableNames);
IBB.Properties.VariableUnits(:) = {'amol'}; %Add units

%CBP80-like cargo
IMs_CBP80 = Yss(height(Yss),:); %Initial masses from steady-state output
IMs_CBP80(strcmp(IMsTable.Properties.VariableNames,'CBP80y')) = Spike*sum(CellVolumeTable{1,:}); %Input mass of spike-in cargo
[~,Y_CBP80] = ode15s(@(t,y)NucCytoShuttleODE(t,y,IMsTable.Properties.VariableNames, ...
    ParameterTable{1,:},ParameterTable.Properties.VariableNames, ...
    CellVolumeTable{1,:},~all(NUPs == 0)),tspan,IMs_CBP80,NucCytoODEopts);
CBP80 = array2table(Y_CBP80,'VariableNames',IMsTable.Properties.VariableNames);
CBP80.Properties.VariableUnits(:) = {'amol'}; %Add units

%NNES-like cargo
IMs_NNES = Yss(height(Yss),:); %Initial masses from steady-state output
IMs_NNES(strcmp(IMsTable.Properties.VariableNames,'NNESy')) = Spike*sum(CellVolumeTable{1,:}); %Input mass of spike-in cargo
[~,Y_NNES] = ode15s(@(t,y)NucCytoShuttleODE(t,y,IMsTable.Properties.VariableNames, ...
    ParameterTable{1,:},ParameterTable.Properties.VariableNames, ...
    CellVolumeTable{1,:},~all(NUPs == 0)),tspan,IMs_NNES,NucCytoODEopts);
NNES = array2table(Y_NNES,'VariableNames',IMsTable.Properties.VariableNames);
NNES.Properties.VariableUnits(:) = {'amol'}; %Add units

%INES-like cargo
IMs_INES = Yss(height(Yss),:); %Initial masses from steady-state output
IMs_INES(strcmp(IMsTable.Properties.VariableNames,'INESy')) = Spike*sum(CellVolumeTable{1,:}); %Input mass of spike-in cargo
[~,Y_INES] = ode15s(@(t,y)NucCytoShuttleODE(t,y,IMsTable.Properties.VariableNames, ...
    ParameterTable{1,:},ParameterTable.Properties.VariableNames, ...
    CellVolumeTable{1,:},~all(NUPs == 0)),tspan,IMs_INES,NucCytoODEopts);
INES = array2table(Y_INES,'VariableNames',IMsTable.Properties.VariableNames);
INES.Properties.VariableUnits(:) = {'amol'}; %Add units

%CNES-like cargo
IMs_CNES = Yss(height(Yss),:); %Initial masses from steady-state output
IMs_CNES(strcmp(IMsTable.Properties.VariableNames,'CNESy')) = Spike*sum(CellVolumeTable{1,:}); %Input mass of spike-in cargo
[~,Y_CNES] = ode15s(@(t,y)NucCytoShuttleODE(t,y,IMsTable.Properties.VariableNames, ...
    ParameterTable{1,:},ParameterTable.Properties.VariableNames, ...
    CellVolumeTable{1,:},~all(NUPs == 0)),tspan,IMs_CNES,NucCytoODEopts);
CNES = array2table(Y_CNES,'VariableNames',IMsTable.Properties.VariableNames);
CNES.Properties.VariableUnits(:) = {'amol'}; %Add units

OutputTable = table(NLS,IBB,CBP80,NNES,INES,CNES,SS);
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