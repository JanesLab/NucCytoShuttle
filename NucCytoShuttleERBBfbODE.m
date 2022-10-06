function dCdt = NucCytoShuttleERBBfbODE(~,I,N,K,Kn,V,NPs,FB)
%NucCytoShuttleODE    Differential equations function for the model
% dCdt = NucCytoShuttleODE(t,I,N,K,Kn,V,NPs)
% t: time (sec); first input, not used in the function
% I: species amounts (amol)
% N: names of species
% K: model parameter (units vary by parameter)
% Kn: names of model parameters
% V: subcellular volumes (pL)
% NPs: logical that evaluates whether nuclear pores are modeled
% FB:  logical that specifies whether ErbB feedback is included

%% Unchanging concentrations
GTP = 470; %µM Table I, Row 11 of Riddick & Macara, J Cell Biol (2005)
GDP = 1.6; %µM Table I, Row 12 of Riddick & Macara, J Cell Biol (2005)

%% Assign an amount variable to each species for ease of interpretation and computational efficiency
Rd = I(strcmp(N,'Rd'));
P = I(strcmp(N,'P'));
RanGap = I(strcmp(N,'RanGap'));
RCC1 = I(strcmp(N,'RCC1'));
Ny = I(strcmp(N,'Ny'));
Ay = I(strcmp(N,'Ay'));
By = I(strcmp(N,'By'));
Sy = I(strcmp(N,'Sy'));
C2y = I(strcmp(N,'C2y'));
Cy = I(strcmp(N,'Cy'));
NRy = I(strcmp(N,'NRy'));
My = I(strcmp(N,'My'));
P3y = I(strcmp(N,'P3y'));
NUP1 = I(strcmp(N,'NUP1'));
NUP2 = I(strcmp(N,'NUP2'));
Rt = I(strcmp(N,'Rt'));
RdNn = I(strcmp(N,'RdNn'));
Rdn = I(strcmp(N,'Rdn'));
Nn = I(strcmp(N,'Nn'));
C2n = I(strcmp(N,'C2n'));
Sn = I(strcmp(N,'Sn'));
Cn = I(strcmp(N,'Cn'));
An = I(strcmp(N,'An'));
Bn = I(strcmp(N,'Bn'));
ABn = I(strcmp(N,'ABn'));
ABCn = I(strcmp(N,'ABCn'));
RtBn = I(strcmp(N,'RtBn'));
RtSAn = I(strcmp(N,'RtSAn'));
C2Bn = I(strcmp(N,'C2Bn'));
ABCFn = I(strcmp(N,'ABCFn'));
ACFn = I(strcmp(N,'ACFn'));
CFn = I(strcmp(N,'CFn'));
C2Fn = I(strcmp(N,'C2Fn'));
C2FBn = I(strcmp(N,'C2FBn'));
RCC1Rd = I(strcmp(N,'RCC1Rd'));
RCC1R = I(strcmp(N,'RCC1R'));
RCC1Rt = I(strcmp(N,'RCC1Rt'));
ABntemp = I(strcmp(N,'ABntemp'));
NRn = I(strcmp(N,'NRn'));
RtNRn = I(strcmp(N,'RtNRn'));
ACn = I(strcmp(N,'ACn'));
RtSn = I(strcmp(N,'RtSn'));
RdNy = I(strcmp(N,'RdNy'));
ABy = I(strcmp(N,'ABy'));
ABCy = I(strcmp(N,'ABCy'));
RtBy = I(strcmp(N,'RtBy'));
RtSAy = I(strcmp(N,'RtSAy'));
C2By = I(strcmp(N,'C2By'));
RtPBy = I(strcmp(N,'RtPBy'));
CFy = I(strcmp(N,'CFy'));
ABCFy = I(strcmp(N,'ABCFy'));
C2Fy = I(strcmp(N,'C2Fy'));
C2FBy = I(strcmp(N,'C2FBy'));
ABCytemp = I(strcmp(N,'ABCytemp'));
ABCFytemp = I(strcmp(N,'ABCFytemp'));
C2Bytemp = I(strcmp(N,'C2Bytemp'));
C2FBytemp = I(strcmp(N,'C2FBytemp'));
RtNRy = I(strcmp(N,'RtNRy'));
Rty = I(strcmp(N,'Rty'));
RtPy = I(strcmp(N,'RtPy'));
ABytemp = I(strcmp(N,'ABytemp'));
RtNRPy = I(strcmp(N,'RtNRPy'));
CBP80y = I(strcmp(N,'CBP80y'));
NESy = I(strcmp(N,'NESy'));
NNESy = I(strcmp(N,'NNESy'));
INESy = I(strcmp(N,'INESy'));
CNESy = I(strcmp(N,'CNESy'));
ABCBP80ytemp = I(strcmp(N,'ABCBP80ytemp'));
ABNNESytemp = I(strcmp(N,'ABNNESytemp'));
ABCNESytemp = I(strcmp(N,'ABCNESytemp'));
ABCBP80y = I(strcmp(N,'ABCBP80y'));
ABCNESy = I(strcmp(N,'ABCNESy'));
BINESy = I(strcmp(N,'BINESy'));
BINESytemp = I(strcmp(N,'BINESytemp'));
ABNNESy = I(strcmp(N,'ABNNESy'));
XPO1y = I(strcmp(N,'XPO1y'));
XPO1Ny = I(strcmp(N,'XPO1Ny'));
XPO1Iy = I(strcmp(N,'XPO1Iy'));
XPO1Cy = I(strcmp(N,'XPO1Cy'));
CBP80n = I(strcmp(N,'CBP80n'));
NESn = I(strcmp(N,'NESn'));
NNESn = I(strcmp(N,'NNESn'));
INESn = I(strcmp(N,'INESn'));
CNESn = I(strcmp(N,'CNESn'));
Mn = I(strcmp(N,'Mn'));
P3n = I(strcmp(N,'P3n'));
MP3n = I(strcmp(N,'MP3n'));
RtMn = I(strcmp(N,'RtMn'));
RtP3n = I(strcmp(N,'RtP3n'));
RtMP3n = I(strcmp(N,'RtMP3n'));
MP3NESn = I(strcmp(N,'MP3NESn'));
MP3NNESn = I(strcmp(N,'MP3NNESn'));
MP3INESn = I(strcmp(N,'MP3INESn'));
MP3CNESn = I(strcmp(N,'MP3CNESn'));
XPO1n = I(strcmp(N,'XPO1n'));
XPO1Cn = I(strcmp(N,'XPO1Cn'));
XPO1In = I(strcmp(N,'XPO1In'));
XPO1Nn = I(strcmp(N,'XPO1Nn'));
ACBP80n = I(strcmp(N,'ACBP80n'));
ACNESn = I(strcmp(N,'ACNESn'));
BINESn = I(strcmp(N,'BINESn'));
ANNESn = I(strcmp(N,'ANNESn'));
ABCBP80n = I(strcmp(N,'ABCBP80n'));
ABCNESn = I(strcmp(N,'ABCNESn'));
ABNNESn = I(strcmp(N,'ABNNESn'));
Byperi = I(strcmp(N,'Byperi'));
Nyperi = I(strcmp(N,'Nyperi'));
Syperi = I(strcmp(N,'Syperi'));
Myperi = I(strcmp(N,'Myperi'));
C2Byperi = I(strcmp(N,'C2Byperi'));
C2FByperi = I(strcmp(N,'C2FByperi'));
RdNyperi = I(strcmp(N,'RdNyperi'));
RtByperi = I(strcmp(N,'RtByperi'));
RtSAyperi = I(strcmp(N,'RtSAyperi'));
AByperi = I(strcmp(N,'AByperi'));
ABCyperi = I(strcmp(N,'ABCyperi'));
ABCFyperi = I(strcmp(N,'ABCFyperi'));
ABCBP80yperi = I(strcmp(N,'ABCBP80yperi'));
ABNNESyperi = I(strcmp(N,'ABNNESyperi'));
ABCNESyperi = I(strcmp(N,'ABCNESyperi'));
BINESyperi = I(strcmp(N,'BINESyperi'));
XPO1Nyperi = I(strcmp(N,'XPO1Nyperi'));
XPO1Iyperi = I(strcmp(N,'XPO1Iyperi'));
XPO1Cyperi = I(strcmp(N,'XPO1Cyperi'));
Bnperi = I(strcmp(N,'Bnperi'));
Nnperi = I(strcmp(N,'Nnperi'));
Snperi = I(strcmp(N,'Snperi'));
Mnperi = I(strcmp(N,'Mnperi'));
C2Bnperi = I(strcmp(N,'C2Bnperi'));
C2FBnperi = I(strcmp(N,'C2FBnperi'));
RdNnperi = I(strcmp(N,'RdNnperi'));
RtBnperi = I(strcmp(N,'RtBnperi'));
RtSAnperi = I(strcmp(N,'RtSAnperi'));
ABnperi = I(strcmp(N,'ABnperi'));
ABCnperi = I(strcmp(N,'ABCnperi'));
ABCFnperi = I(strcmp(N,'ABCFnperi'));
ABCBP80nperi = I(strcmp(N,'ABCBP80nperi'));
ABNNESnperi = I(strcmp(N,'ABNNESnperi'));
ABCNESnperi = I(strcmp(N,'ABCNESnperi'));
BINESnperi = I(strcmp(N,'BINESnperi'));
XPO1Nnperi = I(strcmp(N,'XPO1Nnperi'));
XPO1Inperi = I(strcmp(N,'XPO1Inperi'));
XPO1Cnperi = I(strcmp(N,'XPO1Cnperi'));
NUP1_Bn = I(strcmp(N,'NUP1_Bn'));
NUP1_Nn = I(strcmp(N,'NUP1_Nn'));
NUP1_Sn = I(strcmp(N,'NUP1_Sn'));
NUP1_Mn = I(strcmp(N,'NUP1_Mn'));
NUP1_C2Bn = I(strcmp(N,'NUP1_C2Bn'));
NUP1_C2FBn = I(strcmp(N,'NUP1_C2FBn'));
NUP1_RdNn = I(strcmp(N,'NUP1_RdNn'));
NUP1_RtBn = I(strcmp(N,'NUP1_RtBn'));
NUP1_RtSAn = I(strcmp(N,'NUP1_RtSAn'));
NUP1_ABn = I(strcmp(N,'NUP1_ABn'));
NUP1_ABCn = I(strcmp(N,'NUP1_ABCn'));
NUP1_ABCFn = I(strcmp(N,'NUP1_ABCFn'));
NUP1_ABCBP80n = I(strcmp(N,'NUP1_ABCBP80n'));
NUP1_ABNNESn = I(strcmp(N,'NUP1_ABNNESn'));
NUP1_ABCNESn = I(strcmp(N,'NUP1_ABCNESn'));
NUP1_BINESn = I(strcmp(N,'NUP1_BINESn'));
NUP1_XPO1Nn = I(strcmp(N,'NUP1_XPO1Nn'));
NUP1_XPO1In = I(strcmp(N,'NUP1_XPO1In'));
NUP1_XPO1Cn = I(strcmp(N,'NUP1_XPO1Cn'));
NUP1_By = I(strcmp(N,'NUP1_By'));
NUP1_Ny = I(strcmp(N,'NUP1_Ny'));
NUP1_Sy = I(strcmp(N,'NUP1_Sy'));
NUP1_My = I(strcmp(N,'NUP1_My'));
NUP1_C2By = I(strcmp(N,'NUP1_C2By'));
NUP1_C2FBy = I(strcmp(N,'NUP1_C2FBy'));
NUP1_RdNy = I(strcmp(N,'NUP1_RdNy'));
NUP1_RtBy = I(strcmp(N,'NUP1_RtBy'));
NUP1_RtSAy = I(strcmp(N,'NUP1_RtSAy'));
NUP1_ABy = I(strcmp(N,'NUP1_ABy'));
NUP1_ABCy = I(strcmp(N,'NUP1_ABCy'));
NUP1_ABCFy = I(strcmp(N,'NUP1_ABCFy'));
NUP1_ABCBP80y = I(strcmp(N,'NUP1_ABCBP80y'));
NUP1_ABNNESy = I(strcmp(N,'NUP1_ABNNESy'));
NUP1_ABCNESy = I(strcmp(N,'NUP1_ABCNESy'));
NUP1_BINESy = I(strcmp(N,'NUP1_BINESy'));
NUP1_XPO1Ny = I(strcmp(N,'NUP1_XPO1Ny'));
NUP1_XPO1Iy = I(strcmp(N,'NUP1_XPO1Iy'));
NUP1_XPO1Cy = I(strcmp(N,'NUP1_XPO1Cy'));
ABP3y = I(strcmp(N,'ABP3y'));
ABP3ytemp = I(strcmp(N,'ABP3ytemp'));
ABP3yperi = I(strcmp(N,'ABP3yperi'));
NUP1_ABP3y = I(strcmp(N,'NUP1_ABP3y'));
ABP3n = I(strcmp(N,'ABP3n'));
ABP3nperi = I(strcmp(N,'ABP3nperi'));
NUP1_ABP3n = I(strcmp(N,'NUP1_ABP3n'));
AP3n = I(strcmp(N,'AP3n'));
NUP2_Bn = I(strcmp(N,'NUP2_Bn'));
NUP2_Nn = I(strcmp(N,'NUP2_Nn'));
NUP2_Sn = I(strcmp(N,'NUP2_Sn'));
NUP2_Mn = I(strcmp(N,'NUP2_Mn'));
NUP2_RdNn = I(strcmp(N,'NUP2_RdNn'));
NUP2_RtBn = I(strcmp(N,'NUP2_RtBn'));
NUP2_RtSAn = I(strcmp(N,'NUP2_RtSAn'));
NUP2_ABn = I(strcmp(N,'NUP2_ABn'));
NUP2_By = I(strcmp(N,'NUP2_By'));
NUP2_Ny = I(strcmp(N,'NUP2_Ny'));
NUP2_Sy = I(strcmp(N,'NUP2_Sy'));
NUP2_My = I(strcmp(N,'NUP2_My'));
NUP2_RdNy = I(strcmp(N,'NUP2_RdNy'));
NUP2_RtBy = I(strcmp(N,'NUP2_RtBy'));
NUP2_RtSAy = I(strcmp(N,'NUP2_RtSAy'));
NUP2_ABy = I(strcmp(N,'NUP2_ABy'));
NUP2_ABP3y = I(strcmp(N,'NUP2_ABP3y'));
NUP2_ABP3n = I(strcmp(N,'NUP2_ABP3n'));

%% Assign rate parameters for ease of interpretation and compuational efficiency
P_R = K(strcmp(Kn,'P_R'));
P_N = K(strcmp(Kn,'P_N'));
P_RdN = K(strcmp(Kn,'P_RdN'));
P_B = K(strcmp(Kn,'P_B'));
P_A = K(strcmp(Kn,'P_A'));
P_AB = K(strcmp(Kn,'P_AB'));
P_ABC = K(strcmp(Kn,'P_ABC'));
P_C = K(strcmp(Kn,'P_C'));
P_CF = K(strcmp(Kn,'P_CF'));
P_S = K(strcmp(Kn,'P_S'));
P_NR = K(strcmp(Kn,'P_NR'));
P_RtSA = K(strcmp(Kn,'P_RtSA'));
P_RtB = K(strcmp(Kn,'P_RtB'));
P_RtNR = K(strcmp(Kn,'P_RtNR'));
P_C2B = K(strcmp(Kn,'P_C2B'));
P_M = K(strcmp(Kn,'P_M'));
P_XPO1 = K(strcmp(Kn,'P_XPO1'));
kon_RtB = K(strcmp(Kn,'kon_RtB'));
koff_RtB = K(strcmp(Kn,'koff_RtB'));
kon_RdN = K(strcmp(Kn,'kon_RdN'));
koff_RdN = K(strcmp(Kn,'koff_RdN'));
kon_RCC1Rd = K(strcmp(Kn,'kon_RCC1Rd'));
koff_RCC1Rd = K(strcmp(Kn,'koff_RCC1Rd'));
kon_RCC1RGt = K(strcmp(Kn,'kon_RCC1RGt'));
koff_RCC1RGt = K(strcmp(Kn,'koff_RCC1RGt'));
kon_RCC1RGd = K(strcmp(Kn,'kon_RCC1RGd'));
koff_RCC1RGd = K(strcmp(Kn,'koff_RCC1RGd'));
kon_RCC1Rt = K(strcmp(Kn,'kon_RCC1Rt'));
koff_RCC1Rt = K(strcmp(Kn,'koff_RCC1Rt'));
kon1_AB = K(strcmp(Kn,'kon1_AB'));
koff1_AB = K(strcmp(Kn,'koff1_AB'));
kon2_AB = K(strcmp(Kn,'kon2_AB'));
koff2_AB = K(strcmp(Kn,'koff2_AB'));
kon1_ABC = K(strcmp(Kn,'kon1_ABC'));
koff1_ABC = K(strcmp(Kn,'koff1_ABC'));
kon2_ABC = K(strcmp(Kn,'kon2_ABC'));
koff2_ABC = K(strcmp(Kn,'koff2_ABC'));
kon_ABCRt = K(strcmp(Kn,'kon_ABCRt'));
koff_ABCRt = K(strcmp(Kn,'koff_ABCRt'));
kon_RtS = K(strcmp(Kn,'kon_RtS'));
koff_RtS = K(strcmp(Kn,'koff_RtS'));
kon_RtSAC = K(strcmp(Kn,'kon_RtSAC'));
koff_RtSAC = K(strcmp(Kn,'koff_RtSAC'));
kon1_C2B = K(strcmp(Kn,'kon1_C2B'));
koff1_C2B = K(strcmp(Kn,'koff1_C2B'));
kon2_C2B = K(strcmp(Kn,'kon2_C2B'));
koff2_C2B = K(strcmp(Kn,'koff2_C2B'));
kon_C2BRt = K(strcmp(Kn,'kon_C2BRt'));
koff_C2BRt = K(strcmp(Kn,'koff_C2BRt'));
kon_RtSAP = K(strcmp(Kn,'kon_RtSAP'));
koff_RtSAP = K(strcmp(Kn,'koff_RtSAP'));
kon_RtPB = K(strcmp(Kn,'kon_RtPB'));
koff_RtPB = K(strcmp(Kn,'koff_RtPB'));
kon_RtPBA = K(strcmp(Kn,'kon_RtPBA'));
koff_RtPBA = K(strcmp(Kn,'koff_RtPBA'));
kcat_RtSAPRanGap = K(strcmp(Kn,'kcat_RtSAPRanGap'));
Km_RtSAPRanGap = K(strcmp(Kn,'Km_RtSAPRanGap'));
kcat_RtPBRanGap = K(strcmp(Kn,'kcat_RtPBRanGap'));
Km_RtPBRanGap = K(strcmp(Kn,'Km_RtPBRanGap'));
kcat_RtPRanGap = K(strcmp(Kn,'kcat_RtPRanGap'));
Km_RtPRanGap = K(strcmp(Kn,'Km_RtPRanGap'));
kon_NUPB = K(strcmp(Kn,'kon_NUPB'));
koff_NUPB = K(strcmp(Kn,'koff_NUPB'));
kon_NUPN = K(strcmp(Kn,'kon_NUPN'));
koff_NUPN = K(strcmp(Kn,'koff_NUPN'));
kon_NUPS = K(strcmp(Kn,'kon_NUPS'));
koff_NUPS = K(strcmp(Kn,'koff_NUPS'));
kon_NUPM = K(strcmp(Kn,'kon_NUPM'));
koff_NUPM = K(strcmp(Kn,'koff_NUPM'));
kon_RtP3 = K(strcmp(Kn,'kon_RtP3'));
koff_RtP3 = K(strcmp(Kn,'koff_RtP3'));
kon_MP3 = K(strcmp(Kn,'kon_MP3'));
koff_MP3 = K(strcmp(Kn,'koff_MP3'));
kon_MP3NES = K(strcmp(Kn,'kon_MP3NES'));
koff_MP3NES = K(strcmp(Kn,'koff_MP3NES'));
kon_MP3Rt = K(strcmp(Kn,'kon_MP3Rt'));
koff_MP3Rt = K(strcmp(Kn,'koff_MP3Rt'));
kon_MP3NESRt = K(strcmp(Kn,'kon_MP3NESRt'));
koff_MP3NESRt = K(strcmp(Kn,'koff_MP3NESRt'));
kon_RtMP3NES = K(strcmp(Kn,'kon_RtMP3NES'));
koff_RtMP3NES = K(strcmp(Kn,'koff_RtMP3NES'));
kon_MP3NESRtP = K(strcmp(Kn,'kon_MP3NESRtP'));
koff_MP3NESRtP = K(strcmp(Kn,'koff_MP3NESRtP'));
kon1_ABCBP80 = K(strcmp(Kn,'kon1_ABCBP80'));
koff1_ABCBP80 = K(strcmp(Kn,'koff1_ABCBP80'));
kon2_ABCBP80 = K(strcmp(Kn,'kon2_ABCBP80'));
koff2_ABCBP80 = K(strcmp(Kn,'koff2_ABCBP80'));

%% Assign subcompartment volumes
VolCyto = V(1);
VolNuc = V(2);
VolNPC = V(3);

%% Modular rate processes
%Organized by Jarnac code and then appended with Riddick & Macara, Mol Syst 
%Biol (2007) and finally new elaborations in this model; rate processes are 
%cast in a left (reactant)-to-right (product) direction, with bidirectionality 
%indicated by "<–>".

%Translocation through the NPC
CF_CytoPerm = P_CF * CFn / VolNuc - P_CF * CFy / VolCyto; %CFn <–> CFy
C2F_CytoPerm = P_CF * C2Fn / VolNuc - P_CF * C2Fy / VolCyto; %C2Fn <–> C2Fy
C2_CytoPerm = P_C * C2n / VolNuc - P_C * C2y / VolCyto; %C2n <–> C2y
C_CytoPerm = P_C * Cn / VolNuc - P_C * Cy / VolCyto; %Cn <–> Cy
S_CytoPerm = P_S * Sn / VolNuc - P_S * Sy / VolCyto; %Sn -> Sy; not used in NUP model
RtSA_CytoPerm = P_RtSA * RtSAn / VolNuc - P_RtSA * RtSAy / VolCyto; %RtSAn <–> RtSAy; not used in NUP model
N_CytoPerm = P_N * Nn / VolNuc - P_N * Ny / VolCyto; %Nn <–> Ny; not used in NUP model
Rd_NucPerm = P_R * Rd / VolCyto - P_R * Rdn / VolNuc; %Rd <–> Rdn
Rt_CytoPerm = P_R * Rt / VolNuc - P_R * Rty / VolCyto; %Rt <–> Rty
RdN_NucPerm = P_RdN * RdNy / VolCyto - P_RdN * RdNn / VolNuc; %RdNy <–> RdNn; not used in NUP model
A_CytoPerm = P_A * An / VolNuc - P_A * Ay / VolCyto; %An <–> Ay
B_CytoPerm = P_B * Bn / VolNuc - P_B * By / VolCyto; %Bn <–> By
AB_CytoPerm = P_AB * ABn / VolNuc - P_AB * ABy / VolCyto; %ABn <–> ABy
C2B_NucPerm = P_C2B * C2By / VolCyto - P_C2B * C2Bn / VolNuc; %C2By <–> C2Bn
C2FB_NucPerm = P_C2B * C2FBy / VolCyto - P_C2B * C2FBn / VolNuc; %C2FBy <–> C2FBn
ABC_NucPerm = P_ABC * ABCy / VolCyto - P_ABC * ABCn / VolNuc; %ABCy <–> ABCn
ABCF_NucPerm = P_ABC * ABCFy / VolCyto - P_ABC * ABCFn / VolNuc; %ABCFy <–> ABCFn
RtB_NucPerm = P_RtB * RtBy / VolCyto - P_RtB * RtBn / VolNuc; %RtBy -> RtBn
NR_NucPerm = P_NR * NRy / VolCyto - P_NR * NRn / VolNuc; %NRy <–> NRn
RtNR_CytoPerm = P_RtNR * RtNRn / VolNuc - P_RtNR * RtNRy / VolCyto; %RtNRn <–> RtNRy
%Ran shuttling
Rd_Ny_Binding = kon_RdN * (Rd / VolCyto) * (Ny / VolCyto) ...
    - koff_RdN * RdNy; %Rd + Ny <–> RdNy
%RCC1-mediated GTP exchange
RCC1_Rdn_Binding = kon_RCC1Rd * (Rdn / VolNuc) * (RCC1 / VolNuc) ...
    - koff_RCC1Rd * RCC1Rd / VolNuc; %RCC1 + Rdn <–> RCC1Rd
RCC1_RdNn_Exchange = kon_RCC1Rd * (RdNn / VolNuc) * (RCC1 / VolNuc) ...
    - koff_RCC1Rd * (Nn / VolNuc) * (RCC1Rd / VolNuc) * 1e-6; %RCC1 + RdNn <–> Nn + RCC1Rd; can't use koff_RCC1Rd here but Jarnac code does anyways; rescaling rate parameter
RCC1R_GDP_Unbinding = koff_RCC1RGd * RCC1Rd / VolNuc ...
    - kon_RCC1RGd * (RCC1R / VolNuc) * GDP; %RCC1Rd <–> RCC1R + GDPn
RCC1R_GTP_Binding = kon_RCC1RGt * (RCC1R / VolNuc) * GTP ...
    - koff_RCC1RGt * RCC1Rt / VolNuc; %RCC1R + GTP <–> RCC1Rt
RCC1_Rt_Unbinding = koff_RCC1Rt * RCC1Rt / VolNuc ...
    - kon_RCC1Rt * (RCC1 / VolNuc) * (Rt / VolNuc); %RCC1Rt <–> RCC1 + Rt
%Impα-Impβ-NLS complex assembly in the cytoplasm
Ay_By_Binding = kon1_AB * (Ay / VolCyto) * (By / VolCyto) ...
    - koff1_AB * ABytemp / VolCyto; %Ay + By <–> ABytemp
AyBytemp_Conversion = kon2_AB * ABytemp / VolCyto ...
    - koff2_AB * ABy / VolCyto; %ABytemp <–> ABy
ABy_Cy_Binding = kon1_ABC * (ABy / VolCyto) * (Cy / VolCyto) ...
    - koff1_ABC * ABCytemp / VolCyto; %ABy + Cy <–> ABCytemp
ABCytemp_Conversion = kon2_ABC * ABCytemp / VolCyto ...
    - koff2_ABC * ABCy / VolCyto; %ABCytemp <–> ABCy
%Impα-Impβ-NLS complex disassembly in the nucleus
ABCn_Rt_Exchange = kon_ABCRt * (ABCn / VolNuc) * (Rt / VolNuc) ...
    - koff_ABCRt * (ACn / VolNuc) * (RtBn / VolNuc); %ABCn + Rt <–> ACn + RtBn
Rt_Sn_Binding = kon_RtS * (Rt / VolNuc) * (Sn / VolNuc) ...
    - koff_RtS * RtSn / VolNuc; %Rt + Sn <–> RtSn
RtSn_ACn_Exchange = kon_RtSAC * (RtSn / VolNuc) * (ACn / VolNuc) ...
    - koff_RtSAC * (RtSAn / VolNuc) * (Cn / VolNuc); %RtSn + ACn <–> RtSAn + Cn
Sn_Rt_ACn_Exchange = kon_RtSAC * (Sn / VolNuc) * (Rt / VolNuc) * (ACn / VolNuc) * 1e-6 ...
    - koff_RtSAC * (RtSAn / VolNuc) * (Cn / VolNuc); %Sn + Rt + ACn <–> RtSAn + Cn; can't use kon_RtSAC here but Jarnac code does anyways; rescaling rate parameter
%Impα-Impβ-NLS-Fluor pathway
ABy_CFy_Binding = kon1_ABC * (ABy / VolCyto) * (CFy / VolCyto) ...
    - koff1_ABC * ABCFytemp / VolCyto; %ABy + CFy <–> ABCFytemp
ABCFytemp_Conversion = kon2_ABC * ABCFytemp / VolCyto ...
    - koff2_ABC * ABCFy / VolCyto; %ABCFytemp <–> ABCFy
ABCFn_Rt_Exchange = kon_ABCRt * (ABCFn / VolNuc) * (Rt / VolNuc) ...
    - koff_ABCRt * (ACFn / VolNuc) * (RtBn / VolNuc); %ABCFn + Rt <–> ACFn + RtBn
RtSn_ACFn_Exchange = kon_RtSAC * (RtSn / VolNuc) * (ACFn / VolNuc) ...
    - koff_RtSAC * (RtSAn / VolNuc) * (CFn / VolNuc); %RtSn + ACFn <–> RtSAn + CFn; note that the Jarnac code uses 0.1 x koff_RtSAC for unexplained reasons (assumed to be a typo)
Sn_Rt_ACFn_Exchange = kon_RtSAC * (Sn / VolNuc) * (Rt / VolNuc) * (ACFn / VolNuc) * 1e-6 ...
    - koff_RtSAC * (RtSAn / VolNuc) * (CFn / VolNuc); %Sn + Rt + ACFn <–> RtSAn + CFn; can't use kon_RtSAC here but Jarnac code does anyways; rescaling rate parameter
%Disassembly of complexes in the cytoplasm
RtSAy_P_Disassembly = kon_RtSAP * (RtSAy / VolCyto) * (P / VolCyto) ...
    - koff_RtSAP * (RtPy / VolCyto) * (Ay / VolCyto) * (Sy / VolCyto); %RtSAy + P <–> RtPy + Ay + Sy
RtPy_Disassembly = kcat_RtPRanGap * (RanGap / VolCyto) * (RtPy / VolCyto) ...
    / (Km_RtPRanGap + RtPy / VolCyto);
RtBy_P_Binding = kon_RtPB * (RtBy / VolCyto) * (P / VolCyto) ...
    - koff_RtPB * (RtPBy / VolCyto); %RtBy + P <–> RtPBy
RtPBy_Ay_Exchange = kon_RtPBA * (RtPBy / VolCyto) * (Ay / VolCyto) ...
    - koff_RtPBA * (RtPy / VolCyto) * (ABy / VolCyto); %RtPBy + Ay <–> RtPy + ABy
RtPBy_Disassembly = kcat_RtPBRanGap * (RanGap / VolCyto) * (RtPBy / VolCyto) ...
    / (Km_RtPBRanGap + RtPBy / VolCyto); %RtPBy –> Rd + P + By
%Impβ-Cargo assembly in the cytoplasm
By_C2y_Binding = kon1_C2B * (By / VolCyto) * (C2y / VolCyto) ...
    - koff1_C2B * C2Bytemp / VolCyto; %By + C2y <–> C2Bytemp
C2Bytemp_Conversion = kon2_C2B * C2Bytemp / VolCyto ...
    - koff2_C2B * C2By / VolCyto; %C2Bytemp <–> C2By
%Impβ-Cargo disassembly in the nucleus
C2Bn_Rt_Exchange = kon_C2BRt * (C2Bn / VolNuc) * (Rt / VolNuc) ...
    - koff_C2BRt * (C2n / VolNuc) * (RtBn / VolNuc); %C2Bn + Rt <–> C2n + RtBn
%Impβ-IBB-Fluor pathway
By_C2Fy_Binding = kon1_C2B * (By / VolCyto) * (C2Fy / VolCyto) ...
    - koff1_C2B * C2FBytemp / VolCyto; %By + C2Fy <–> C2FBytemp
C2FBytemp_Conversion = kon2_C2B * C2FBytemp / VolCyto ...
    - koff2_C2B * C2FBy / VolCyto; %C2FBytemp <–> C2FBy
C2FBn_Rt_Exchange = kon_C2BRt * (C2FBn / VolNuc) * (Rt / VolNuc) ...
    - koff_C2BRt * (C2Fn / VolNuc) * (RtBn / VolNuc); %C2FBn + Rt <–> C2Fn + RtBn
%Carrier cycling
An_Bn_Binding = kon1_AB * (An / VolNuc) * (Bn / VolNuc) ...
    - koff1_AB * ABntemp / VolNuc; %An + Bn <–> ABntemp
ABntemp_Conversion = kon2_AB * ABntemp / VolNuc ...
    - koff2_AB * ABn / VolNuc; %ABntemp <–> ABn
Rt_Bn_Binding = kon_RtB * (Rt / VolNuc) * (Bn / VolNuc) ...
    - koff_RtB * RtBn / VolNuc; %Rt + Bn <–> RtBn
ABn_Rt_Exchange = kon_C2BRt * (ABn / VolNuc) * (Rt / VolNuc) ...
    - koff_C2BRt * (An / VolNuc) * (RtBn / VolNuc); %ABn + Rt <–> An + RtBn; taking Impα as a C2B-like cargo
%Cycling of other transport receptors
Rt_NRn_Binding = kon_RtB * (Rt / VolNuc) * (NRn / VolNuc) ...
    - koff_RtB * RtNRn / VolNuc; %Rt + NRn <–> RtNRn; other transport receptors assumed to bind Rt like Impβ
RtNRy_P_Binding = kon_RtPB * (RtNRy / VolCyto) * (P / VolCyto) ...
    - koff_RtPB * RtNRPy / VolCyto; %RtNRy + P <–> RtNRPy; other transport receptors assumed to bind P like Impβ
RtNRPy_Disassembly = kcat_RtPRanGap * (RanGap / VolCyto) * (RtNRPy / VolCyto) ...
    / (Km_RtPRanGap + RtNRPy / VolCyto); %RtNRPy –> Rd + P + NRy
%Reactions due to RanGTP leakage into the cytoplasm
Rt_By_Binding = kon_RtB * (Rty / VolCyto) * (By / VolCyto) ...
    - koff_RtB * RtBy / VolCyto; %Rty + By <–> RtBy
Sy_Rty_Ay_Binding = kon_RtSAC * (Sy / VolCyto) * (Rty / VolCyto) * (Ay / VolCyto) * 1e-6 ...
    - koff_RtSAC * RtSAy / VolCyto; %Sy + Rty + Ay <–> RtSAy; can't use kon_RtSAC here but Jarnac code does anyways; rescaling rate parameter
Rt_Hydrolysis = kcat_RtPRanGap * (RanGap / VolCyto) * (Rty / VolCyto) ...
    / (Km_RtPRanGap + Rty / VolCyto); %Rty –> Rd; free Rty assumed to be hydrolyzed by RanGap like RtP
Rt_P_Binding = kon_RtPB * (Rty / VolCyto) * (P / VolCyto) ...
    - koff_RtPB * RtPy / VolCyto; %Rty + P <–> RtPy; free Rty assumed to bind P like RtB
%Reactions involving RanBP3 (tabulated)
Rt_P3n_Binding = kon_RtP3 * (Rt / VolNuc) * (P3n / VolNuc) ...
    - koff_RtP3 * RtP3n / VolNuc; %Rt + P3n <–> RtP3n
Mn_P3n_Binding = kon_MP3 * (Mn / VolNuc) * (P3n / VolNuc) ...
    - koff_MP3 * MP3n / VolNuc; %Mn + P3n <–> MP3n
MP3n_NESn_Binding = kon_MP3NES * (MP3n / VolNuc) * (NESn / VolNuc) ...
    - koff_MP3NES * MP3NESn / VolNuc; %MP3n + NESn <–> MP3NESn
Rt_MP3n_Binding = kon_MP3Rt * (Rt / VolNuc) * (MP3n / VolNuc) ...
    - koff_MP3Rt * RtMP3n / VolNuc; %Rt + MP3n <–> RtMP3n
Rt_MP3NESn_Binding = kon_MP3NESRt * (Rt / VolNuc) * (MP3NESn / VolNuc) ...
    - koff_MP3NESRt * XPO1n / VolNuc; %Rt + MP3NESn <–> XPO1n; XPO1 is shorthand for RtMP3NESn
RtMP3n_NESn_Binding = kon_RtMP3NES * (RtMP3n / VolNuc) * (NESn / VolNuc) ...
    - koff_RtMP3NES * XPO1n / VolNuc; %RtMP3n + NESn <–> XPO1n; XPO1 is shorthand for RtMP3NESn
XPO1y_P_Disassembly = kon_MP3NESRtP * (XPO1y / VolCyto) * (P / VolCyto) ...
    - koff_MP3NESRtP * (RtPy / VolCyto) * (My / VolCyto) ...
    * (P3y / VolCyto) * (NESy / VolCyto); %XPO1y + P <–> RtPy + My + P3y + NESy
%Reactions involving RanBP3 (assumed)
Rt_Mn_Binding = kon_MP3 * (Rt / VolNuc) * (Mn / VolNuc) ...
    - koff_MP3 * RtMn / VolNuc; %Rt + Mn <–> RtMn; assumes same binding kinetics as Mn and P3n
RtMn_P3n_Binding = kon_MP3 * (RtMn / VolNuc) * (P3n / VolNuc) ...
    - koff_MP3 * RtMP3n / VolNuc; %RtMn + P3n <–> RtMP3n; assumes same binding kinetics as Mn and P3n
RtP3n_Mn_Binding = kon_MP3 * (RtP3n / VolNuc) * (Mn / VolNuc) ...
    - koff_MP3 * RtMP3n / VolNuc; %RtP3n + Mn <–> RtMP3n; assumes same binding kinetics as Mn and P3n
%Reactions involving NES cargo
NES_CytoPerm = P_CF * NESn / VolNuc - P_CF * NESy / VolCyto; %NESn <–> NESy
XPO1_CytoPerm = P_XPO1 * XPO1n / VolNuc - P_XPO1 * XPO1y / VolCyto; %XPO1n <–> XPO1y
%Reactions involving CBP80
CBP80_CytoPerm = P_CF * CBP80n / VolNuc - P_CF * CBP80y / VolCyto; %CBP80n <–> CBP80y
ABCBP80_NucPerm = P_ABC * ABCBP80y / VolCyto - P_ABC * ABCBP80n / VolNuc; %ABCBP80y <–> ABCBP80n
RtSn_ACBP80n_Exchange = kon_RtSAC * (RtSn / VolNuc) * (ACBP80n / VolNuc) ...
    - koff_RtSAC * (RtSAn / VolNuc) * (CBP80n / VolNuc); %RtSn + ACBP80n <–> RtSAn + CBP80n
Sn_Rt_ACBP80n_Exchange = kon_RtSAC * (Sn / VolNuc) * (Rt / VolNuc) * (ACBP80n / VolNuc) * 1e-6 ...
    - koff_RtSAC * (RtSAn / VolNuc) * (CBP80n / VolNuc); %Sn + Rt + ACBP80n <–> RtSAn + CBP80n; can't use kon_RtSAC here but Jarnac code does anyways; rescaling rate parameter
ABCBP80n_Rt_Exchange = kon_ABCRt * (ABCBP80n / VolNuc) * (Rt / VolNuc) ...
    - koff_ABCRt * (ACBP80n / VolNuc) * (RtBn / VolNuc); %ABCBP80n + Rt <–> ACBP80n + RtBn
ABy_CBP80y_Binding = kon1_ABCBP80 * (ABy / VolCyto) * (CBP80y / VolCyto) ...
    - koff1_ABCBP80 * ABCBP80ytemp / VolCyto; %ABy + CBP80y <–> ABCBP80ytemp
ABCBP80ytemp_Conversion = kon2_ABCBP80 * ABCBP80ytemp / VolCyto ...
    - koff2_ABCBP80 * ABCBP80y / VolCyto; %ABCBP80ytemp <–> ABCBP80y
%Reactions involving NNES
NNES_CytoPerm = P_CF * NNESn / VolNuc - P_CF * NNESy / VolCyto; %NNESn <–> NNESy
RtSn_ANNESn_Exchange = kon_RtSAC * (RtSn / VolNuc) * (ANNESn / VolNuc) ...
    - koff_RtSAC * (RtSAn / VolNuc) * (NNESn / VolNuc); %RtSn + ANNESn <–> RtSAn + NNESn
Sn_Rt_ANNESn_Exchange = kon_RtSAC * (Sn / VolNuc) * (Rt / VolNuc) * (ANNESn / VolNuc) * 1e-6 ...
    - koff_RtSAC * (RtSAn / VolNuc) * (NNESn / VolNuc); %Sn + Rt + ANNESn <–> RtSAn + NNESn; can't use kon_RtSAC here but Jarnac code does anyways; rescaling rate parameter
ABNNESn_Rt_Exchange = kon_ABCRt * (ABNNESn / VolNuc) * (Rt / VolNuc) ...
    - koff_ABCRt * (ANNESn / VolNuc) * (RtBn / VolNuc); %ABNNESn + Rt <–> ANNESn + RtBn
ABy_NNESy_Binding = kon1_ABC * (ABy / VolCyto) * (NNESy / VolCyto) ...
    - koff1_ABC * ABNNESytemp / VolCyto; %ABy + NNESy <–> ABNNESytemp
ABNNESytemp_Conversion = kon2_ABC * ABNNESytemp / VolCyto ...
    - koff2_ABC * ABNNESy / VolCyto; %ABNNESytemp <–> ABNNESy
MP3n_NNESn_Binding = kon_MP3NES * (MP3n / VolNuc) * (NNESn / VolNuc) ...
    - koff_MP3NES * MP3NNESn / VolNuc; %MP3n + NNESn <–> MP3NNESn
Rt_MP3NNESn_Binding = kon_MP3NESRt * (Rt / VolNuc) * (MP3NNESn / VolNuc) ...
    - koff_MP3NESRt * XPO1Nn / VolNuc; %Rt + MP3NNESn <–> XPO1Nn; XPO1Nn is shorthand for RtMP3NNESn
RtMP3n_NNESn_Binding = kon_RtMP3NES * (RtMP3n / VolNuc) * (NNESn / VolNuc) ...
    - koff_RtMP3NES * XPO1Nn / VolNuc; %RtMP3n + NNESn <–> XPO1Nn; XPO1Nn is shorthand for RtMP3NNESn
XPO1Ny_P_Disassembly = kon_MP3NESRtP * (XPO1Ny / VolCyto) * (P / VolCyto) ...
    - koff_MP3NESRtP * (RtPy / VolCyto) * (My / VolCyto) ...
    * (P3y / VolCyto) * (NNESy / VolCyto); %XPO1Ny + P <–> RtPy + My + P3y + NNESy
%Reactions involving CNES
CNES_CytoPerm = P_CF * CNESn / VolNuc - P_CF * CNESy / VolCyto; %CNESn <–> CNESy
RtSn_ACNESn_Exchange = kon_RtSAC * (RtSn / VolNuc) * (ACNESn / VolNuc) ...
    - koff_RtSAC * (RtSAn / VolNuc) * (CNESn / VolNuc); %RtSn + ACNESn <–> RtSAn + CNESn
Sn_Rt_ACNESn_Exchange = kon_RtSAC * (Sn / VolNuc) * (Rt / VolNuc) * (ACNESn / VolNuc) * 1e-6 ...
    - koff_RtSAC * (RtSAn / VolNuc) * (CNESn / VolNuc); %Sn + Rt + ACNESn <–> RtSAn + CNESn; can't use kon_RtSAC here but Jarnac code does anyways; rescaling rate parameter
ABCNESn_Rt_Exchange = kon_ABCRt * (ABCNESn / VolNuc) * (Rt / VolNuc) ...
    - koff_ABCRt * (ACNESn / VolNuc) * (RtBn / VolNuc); %ABCNESn + Rt <–> ACNESn + RtBn
ABy_CNESy_Binding = kon1_ABCBP80 * (ABy / VolCyto) * (CNESy / VolCyto) ...
    - koff1_ABCBP80 * ABCNESytemp / VolCyto; %ABy + CNESy <–> ABCNESytemp
ABCNESytemp_Conversion = kon2_ABCBP80 * ABCNESytemp / VolCyto ...
    - koff2_ABCBP80 * ABCNESy / VolCyto; %ABCNESytemp <–> ABCNESy
MP3n_CNESn_Binding = kon_MP3NES * (MP3n / VolNuc) * (CNESn / VolNuc) ...
    - koff_MP3NES * MP3CNESn / VolNuc; %MP3n + CNESn <–> MP3CNESn
Rt_MP3CNESn_Binding = kon_MP3NESRt * (Rt / VolNuc) * (MP3CNESn / VolNuc) ...
    - koff_MP3NESRt * XPO1Cn / VolNuc; %Rt + MP3CNESn <–> XPO1Cn; XPO1Cn is shorthand for RtMP3CNESn
RtMP3n_CNESn_Binding = kon_RtMP3NES * (RtMP3n / VolNuc) * (CNESn / VolNuc) ...
    - koff_RtMP3NES * XPO1Cn / VolNuc; %RtMP3n + CNESn <–> XPO1Cn; XPO1Cn is shorthand for RtMP3CNESn
XPO1Cy_P_Disassembly = kon_MP3NESRtP * (XPO1Cy / VolCyto) * (P / VolCyto) ...
    - koff_MP3NESRtP * (RtPy / VolCyto) * (My / VolCyto) ...
    * (P3y / VolCyto) * (CNESy / VolCyto); %XPO1Cy + P <–> RtPy + My + P3y + CNESy
%Reactions involving INES
INES_CytoPerm = P_CF * INESn / VolNuc - P_CF * INESy / VolCyto; %INESn <–> INESy
BINESn_Rt_Exchange = kon_C2BRt * (BINESn / VolNuc) * (Rt / VolNuc) ...
    - koff_C2BRt * (INESn / VolNuc) * (RtBn / VolNuc); %BINESn + Rt <–> INESn + RtBn
By_INESy_Binding = kon1_C2B * (By / VolCyto) * (INESy / VolCyto) ...
    - koff1_C2B * BINESytemp / VolCyto; %By + INESy <–> BINESytemp
BINESy_Conversion = kon2_C2B * BINESytemp / VolCyto ...
    - koff2_C2B * BINESy / VolCyto; %BINESytemp <–> BINESy
MP3n_INESn_Binding = kon_MP3NES * (MP3n / VolNuc) * (INESn / VolNuc) ...
    - koff_MP3NES * MP3INESn / VolNuc; %MP3n + INESn <–> MP3INESn
Rt_MP3INESn_Binding = kon_MP3NESRt * (Rt / VolNuc) * (MP3INESn / VolNuc) ...
    - koff_MP3NESRt * XPO1In / VolNuc; %Rt + MP3INESn <–> XPO1In; XPO1I is shorthand for RtMP3INESn
RtMP3n_INESn_Binding = kon_RtMP3NES * (RtMP3n / VolNuc) * (INESn / VolNuc) ...
    - koff_RtMP3NES * XPO1In / VolNuc; %RtMP3n + INESn <–> XPO1In; XPO1I is shorthand for RtMP3INESn
XPO1Iy_P_Disassembly = kon_MP3NESRtP * (XPO1Iy / VolCyto) * (P / VolCyto) ...
    - koff_MP3NESRtP * (RtPy / VolCyto) * (My / VolCyto) ...
    * (P3y / VolCyto) * (INESy / VolCyto); %XPO1Iy + P <–> RtPy + My + P3y + INESy
%Reactions involving cytoplasmic transport to the perinucleus
By_PeriPerm = P_B * By / VolCyto - P_B * Byperi / VolNPC; %By <–> Byperi
Ny_PeriPerm = P_N * Ny / VolCyto - P_N * Nyperi / VolNPC; %Ny <–> Nyperi
Sy_PeriPerm = P_S * Sy / VolCyto - P_S * Syperi / VolNPC; %Sy <–> Syperi
My_PeriPerm = P_M * My / VolCyto - P_M * Myperi / VolNPC; %My <–> Myperi
RtSAy_PeriPerm = P_RtSA * RtSAy / VolCyto - P_RtSA * RtSAyperi / VolNPC; %RtSAy <–> RtSAyperi
ABy_PeriPerm = P_AB * ABy / VolCyto - P_AB * AByperi / VolNPC; %ABy <–> AByperi
C2By_PeriPerm = P_C2B * C2By / VolCyto - P_C2B * C2Byperi / VolNPC; %C2By <–> C2Byperi
C2FBy_PeriPerm = P_C2B * C2FBy / VolCyto - P_C2B * C2FByperi / VolNPC; %C2FBy <–> C2FByperi
ABCy_PeriPerm = P_ABC * ABCy / VolCyto - P_ABC * ABCyperi / VolNPC; %ABCy <–> ABCyperi
ABCFy_PeriPerm = P_ABC * ABCFy / VolCyto - P_ABC * ABCFyperi / VolNPC; %ABCFy <–> ABCFyperi
RdNy_PeriPerm = P_RdN * RdNy / VolCyto - P_RdN * RdNyperi / VolNPC; %RdNy <–> RdNyperi
RtBy_PeriPerm = P_RtB * RtBy / VolCyto - P_RtB * RtByperi / VolNPC; %RtBy <–> RtByperi
ABCBP80y_PeriPerm = P_ABC * ABCBP80y / VolCyto - P_ABC * ABCBP80yperi / VolNPC; %ABCBP80y <–> ABCBP80yperi
XPO1Ny_PeriPerm = P_XPO1 * XPO1Ny / VolCyto - P_XPO1 * XPO1Nyperi / VolNPC; %XPO1Ny <–> XPO1Nyperi
ABNNESy_PeriPerm = P_ABC * ABNNESy / VolCyto - P_ABC * ABNNESyperi / VolNPC; %ABNNESy <–> ABNNESyperi
XPO1Cy_PeriPerm = P_XPO1 * XPO1Cy / VolCyto - P_XPO1 * XPO1Cyperi / VolNPC; %XPO1Cy <–> XPO1Cyperi
ABCNESy_PeriPerm = P_ABC * ABCNESy / VolCyto - P_ABC * ABCNESyperi / VolNPC; %ABCNESy <–> ABCNESyperi
XPO1Iy_PeriPerm = P_XPO1 * XPO1Iy / VolCyto - P_XPO1 * XPO1Iyperi / VolNPC; %XPO1Iy <–> XPO1Iyperi
BINESy_PeriPerm = P_C2B * BINESy / VolCyto - P_C2B * BINESyperi / VolNPC; %BINESy <–> BINESyperi
%Reactions involving perinuclear binding to the primary NUP (cytoplasmic side)
Byperi_NUP1_Binding = kon_NUPB * (Byperi / VolNPC) * (NUP1 / VolNPC); %Byperi + NUP1 –>  NUP1_By
Nyperi_NUP1_Binding = kon_NUPN * (Nyperi / VolNPC) * (NUP1 / VolNPC); %Nyperi + NUP1 –>  NUP1_Ny
Syperi_NUP1_Binding  = kon_NUPS * (Syperi / VolNPC) * (NUP1 / VolNPC); %Syperi + NUP1 –>  NUP1_Sy
Myperi_NUP1_Binding = kon_NUPM * (Myperi / VolNPC) * (NUP1 / VolNPC); %Myperi + NUP1 –>  NUP1_My
RtSAyperi_NUP1_Binding = kon_NUPS * (RtSAyperi / VolNPC) * (NUP1 / VolNPC); %RtSAyperi + NUP1 –>  NUP1_RtSAy
AByperi_NUP1_Binding = kon_NUPB * (AByperi / VolNPC) * (NUP1 / VolNPC); %AByperi + NUP1 –>  NUP1_ABy
C2Byperi_NUP1_Binding = kon_NUPB * (C2Byperi / VolNPC) * (NUP1 / VolNPC); %C2Byperi + NUP1 –>  NUP1_C2By
C2FByperi_NUP1_Binding = kon_NUPB * (C2FByperi / VolNPC) * (NUP1 / VolNPC); %C2FByperi + NUP1 –>  NUP1_C2FBy
ABCyperi_NUP1_Binding = kon_NUPB * (ABCyperi / VolNPC) * (NUP1 / VolNPC); %ABCyperi + NUP1 –>  NUP1_ABCy
ABCFyperi_NUP1_Binding = kon_NUPB * (ABCFyperi / VolNPC) * (NUP1 / VolNPC); %ABCFyperi + NUP1 –>  NUP1_ABCFy
RdNyperi_NUP1_Binding = kon_NUPN * (RdNyperi / VolNPC) * (NUP1 / VolNPC); %RdNyperi + NUP1 –>  NUP1_RdNy
RtByperi_NUP1_Binding = kon_NUPB * (RtByperi / VolNPC) * (NUP1 / VolNPC); %RtByperi + NUP1 –>  NUP1_RtBy
ABCBP80yperi_NUP1_Binding = kon_NUPB * (ABCBP80yperi / VolNPC) * (NUP1 / VolNPC); %ABCBP80yperi + NUP1 –>  NUP1_ABCBP80y
XPO1Nyperi_NUP1_Binding = kon_NUPB * (XPO1Nyperi / VolNPC) * (NUP1 / VolNPC); %XPO1Nyperi + NUP1 –>  NUP1_XPO1Ny
ABNNESyperi_NUP1_Binding = kon_NUPB * (ABNNESyperi / VolNPC) * (NUP1 / VolNPC); %ABNNESyperi + NUP1 –>  NUP1_ABNNESy
XPO1Cyperi_NUP1_Binding = kon_NUPB * (XPO1Cyperi / VolNPC) * (NUP1 / VolNPC); %XPO1Cyperi + NUP1 –>  NUP1_XPO1Cy
ABCNESyperi_NUP1_Binding = kon_NUPB * (ABCNESyperi / VolNPC) * (NUP1 / VolNPC); %ABCNESyperi + NUP1 –>  NUP1_ABCNESy
XPO1Iyperi_NUP1_Binding = kon_NUPB * (XPO1Iyperi / VolNPC) * (NUP1 / VolNPC); %XPO1Iyperi + NUP1 –>  NUP1_XPO1Iy
BINESyperi_NUP1_Binding = kon_NUPB * (BINESyperi / VolNPC) * (NUP1 / VolNPC); %BINESyperi + NUP1 –>  NUP1_BINESy
%Reactions involving primary NUP release into the nucleus
NUP1_By_Dissociation = koff_NUPB * NUP1_By / VolNPC; %NUP1_By –> NUP1 + Bn
NUP1_Ny_Dissociation = koff_NUPN * NUP1_Ny / VolNPC; %NUP1_Ny –> NUP1 + Nn
NUP1_Sy_Dissociation = koff_NUPS * NUP1_Sy / VolNPC; %NUP1_Sy –> NUP1 + Sn
NUP1_My_Dissociation = koff_NUPM * NUP1_My / VolNPC; %NUP1_My –> NUP1 + Mn
NUP1_RtSAy_Dissociation = koff_NUPS * NUP1_RtSAy / VolNPC; %NUP1_RtSAy –> NUP1 + RtSAn
NUP1_ABy_Dissociation = koff_NUPB * NUP1_ABy / VolNPC; %NUP1_ABy –> NUP1 + ABn
NUP1_C2By_Dissociation = koff_NUPB * NUP1_C2By / VolNPC; %NUP1_C2By –> NUP1 + C2Bn
NUP1_C2FBy_Dissociation = koff_NUPB * NUP1_C2FBy / VolNPC; %NUP1_C2FBy –> NUP1 + C2FBn
NUP1_ABCy_Dissociation = koff_NUPB * NUP1_ABCy / VolNPC; %NUP1_ABCy –> NUP1 + ABCn
NUP1_ABCFy_Dissociation = koff_NUPB * NUP1_ABCFy / VolNPC; %NUP1_ABCFy –> NUP1 + ABCFn
NUP1_RdNy_Dissociation = koff_NUPN * NUP1_RdNy / VolNPC; %NUP1_RdNy –> NUP1 + RdNn
NUP1_RtBy_Dissociation = koff_NUPB * NUP1_RtBy / VolNPC; %NUP1_RtBy –> NUP1 + RtBn
NUP1_ABCBP80y_Dissociation = koff_NUPB * NUP1_ABCBP80y / VolNPC; %NUP1_ABCBP80y –> NUP1 + ABCBP80n
NUP1_XPO1Ny_Dissociation = koff_NUPB * NUP1_XPO1Ny / VolNPC; %NUP1_XPO1Ny –> NUP1 + XPO1Nn
NUP1_ABNNESy_Dissociation = koff_NUPB * NUP1_ABNNESy / VolNPC; %NUP1_ABNNESy –> NUP1 + ABNNESn
NUP1_XPO1Cy_Dissociation = koff_NUPB * NUP1_XPO1Cy / VolNPC; %NUP1_XPO1Cy –> NUP1 + XPO1Cn
NUP1_ABCNESy_Dissociation = koff_NUPB * NUP1_ABCNESy / VolNPC; %NUP1_ABCNESy –> NUP1 + ABCNESn
NUP1_XPO1Iy_Dissociation = koff_NUPB * NUP1_XPO1Iy / VolNPC; %NUP1_XPO1Iy –> NUP1 + XPO1In
NUP1_BINESy_Dissociation = koff_NUPB * NUP1_BINESy / VolNPC; %NUP1_BINESy –> NUP1 + BINESn
%Reactions involving nuclear transport to the perinucleus
Bn_PeriPerm = P_B * Bn / VolNuc - P_B * Bnperi / VolNPC; %Bn <–> Bnperi
Nn_PeriPerm = P_N * Nn / VolNuc - P_N * Nnperi / VolNPC; %Nn <–> Nnperi
Sn_PeriPerm = P_S * Sn / VolNuc - P_S * Snperi / VolNPC; %Sn <–> Snperi
Mn_PeriPerm = P_M * Mn / VolNuc - P_M * Mnperi / VolNPC; %Mn <–> Mnperi
RtSAn_PeriPerm = P_RtSA * RtSAn / VolNuc - P_RtSA * RtSAnperi / VolNPC; %RtSAn <–> RtSAnperi
ABn_PeriPerm = P_AB * ABn / VolNuc - P_AB * ABnperi / VolNPC; %ABn <–> ABnperi
C2Bn_PeriPerm = P_C2B * C2Bn / VolNuc - P_C2B * C2Bnperi / VolNPC; %C2Bn <–> C2Bnperi
C2FBn_PeriPerm = P_C2B * C2FBn / VolNuc - P_C2B * C2FBnperi / VolNPC; %C2FBn <–> C2FBnperi
ABCn_PeriPerm = P_ABC * ABCn / VolNuc - P_ABC * ABCnperi / VolNPC; %ABCn <–> ABCnperi
ABCFn_PeriPerm = P_ABC * ABCFn / VolNuc - P_ABC * ABCFnperi / VolNPC; %ABCFn <–> ABCFnperi
RdNn_PeriPerm = P_RdN * RdNn / VolNuc - P_RdN * RdNnperi / VolNPC; %RdNn <–> RdNnperi
RtBn_PeriPerm = P_RtB * RtBn / VolNuc - P_RtB * RtBnperi / VolNPC; %RtBn <–> RtBnperi
ABCBP80n_PeriPerm = P_ABC * ABCBP80n / VolNuc - P_ABC * ABCBP80nperi / VolNPC; %ABCBP80n <–> ABCBP80nperi
XPO1Nn_PeriPerm = P_XPO1 * XPO1Nn / VolNuc - P_XPO1 * XPO1Nnperi / VolNPC; %XPO1Nn <–> XPO1Nnperi
ABNNESn_PeriPerm = P_ABC * ABNNESn / VolNuc - P_ABC * ABNNESnperi / VolNPC; %ABNNESn <–> ABNNESnperi
XPO1Cn_PeriPerm = P_XPO1 * XPO1Cn / VolNuc - P_XPO1 * XPO1Cnperi / VolNPC; %XPO1Cn <–> XPO1Cnperi
ABCNESn_PeriPerm = P_ABC * ABCNESn / VolNuc - P_ABC * ABCNESnperi / VolNPC; %ABCNESn <–> ABCNESnperi
XPO1In_PeriPerm = P_XPO1 * XPO1In / VolNuc - P_XPO1 * XPO1Inperi / VolNPC; %XPO1In <–> XPO1Inperi
BINESn_PeriPerm = P_C2B * BINESn / VolNuc - P_C2B * BINESnperi / VolNPC; %BINESn <–> BINESnperi
%Reactions involving perinuclear binding to the primary NUP (nuclear side)
Bnperi_NUP1_Binding = kon_NUPB * (Bnperi / VolNPC) * (NUP1 / VolNPC); %Bnperi + NUP1 –>  NUP1_Bn
Nnperi_NUP1_Binding = kon_NUPN * (Nnperi / VolNPC) * (NUP1 / VolNPC); %Nnperi + NUP1 –>  NUP1_Nn
Snperi_NUP1_Binding  = kon_NUPS * (Snperi / VolNPC) * (NUP1 / VolNPC); %Snperi + NUP1 –>  NUP1_Sn
Mnperi_NUP1_Binding = kon_NUPM * (Mnperi / VolNPC) * (NUP1 / VolNPC); %Mnperi + NUP1 –>  NUP1_Mn
RtSAnperi_NUP1_Binding = kon_NUPS * (RtSAnperi / VolNPC) * (NUP1 / VolNPC); %RtSAnperi + NUP1 –>  NUP1_RtSAn
ABnperi_NUP1_Binding = kon_NUPB * (ABnperi / VolNPC) * (NUP1 / VolNPC); %ABnperi + NUP1 –>  NUP1_ABn
C2Bnperi_NUP1_Binding = kon_NUPB * (C2Bnperi / VolNPC) * (NUP1 / VolNPC); %C2Bnperi + NUP1 –>  NUP1_C2Bn
C2FBnperi_NUP1_Binding = kon_NUPB * (C2FBnperi / VolNPC) * (NUP1 / VolNPC); %C2FBnperi + NUP1 –>  NUP1_C2FBn
ABCnperi_NUP1_Binding = kon_NUPB * (ABCnperi / VolNPC) * (NUP1 / VolNPC); %ABCnperi + NUP1 –>  NUP1_ABCn
ABCFnperi_NUP1_Binding = kon_NUPB * (ABCFnperi / VolNPC) * (NUP1 / VolNPC); %ABCFnperi + NUP1 –>  NUP1_ABCFn
RdNnperi_NUP1_Binding = kon_NUPN * (RdNnperi / VolNPC) * (NUP1 / VolNPC); %RdNnperi + NUP1 –>  NUP1_RdNn
RtBnperi_NUP1_Binding = kon_NUPB * (RtBnperi / VolNPC) * (NUP1 / VolNPC); %RtBnperi + NUP1 –>  NUP1_RtBn
ABCBP80nperi_NUP1_Binding = kon_NUPB * (ABCBP80nperi / VolNPC) * (NUP1 / VolNPC); %ABCBP80nperi + NUP1 –>  NUP1_ABCBP80n
XPO1Nnperi_NUP1_Binding = kon_NUPB * (XPO1Nnperi / VolNPC) * (NUP1 / VolNPC); %XPO1Nnperi + NUP1 –>  NUP1_XPO1Nn
ABNNESnperi_NUP1_Binding = kon_NUPB * (ABNNESnperi / VolNPC) * (NUP1 / VolNPC); %ABNNESnperi + NUP1 –>  NUP1_ABNNESn
XPO1Cnperi_NUP1_Binding = kon_NUPB * (XPO1Cnperi / VolNPC) * (NUP1 / VolNPC); %XPO1Cnperi + NUP1 –>  NUP1_XPO1Cn
ABCNESnperi_NUP1_Binding = kon_NUPB * (ABCNESnperi / VolNPC) * (NUP1 / VolNPC); %ABCNESnperi + NUP1 –>  NUP1_ABCNESn
XPO1Inperi_NUP1_Binding = kon_NUPB * (XPO1Inperi / VolNPC) * (NUP1 / VolNPC); %XPO1Inperi + NUP1 –>  NUP1_XPO1In
BINESnperi_NUP1_Binding = kon_NUPB * (BINESnperi / VolNPC) * (NUP1 / VolNPC); %BINESnperi + NUP1 –>  NUP1_BINESn
%Reactions involving primary NUP release into the cytoplasm
NUP1_Bn_Dissociation = koff_NUPB * NUP1_Bn / VolNPC; %NUP1_Bn –> NUP1 + By
NUP1_Nn_Dissociation = koff_NUPN * NUP1_Nn / VolNPC; %NUP1_Nn –> NUP1 + Ny
NUP1_Sn_Dissociation = koff_NUPS * NUP1_Sn / VolNPC; %NUP1_Sn –> NUP1 + Sy
NUP1_Mn_Dissociation = koff_NUPM * NUP1_Mn / VolNPC; %NUP1_Mn –> NUP1 + My
NUP1_RtSAn_Dissociation = koff_NUPS * NUP1_RtSAn / VolNPC; %NUP1_RtSAn –> NUP1 + RtSAy
NUP1_ABn_Dissociation = koff_NUPB * NUP1_ABn / VolNPC; %NUP1_ABn –> NUP1 + ABy
NUP1_C2Bn_Dissociation = koff_NUPB * NUP1_C2Bn / VolNPC; %NUP1_C2Bn –> NUP1 + C2By
NUP1_C2FBn_Dissociation = koff_NUPB * NUP1_C2FBn / VolNPC; %NUP1_C2FBn –> NUP1 + C2FBy
NUP1_ABCn_Dissociation = koff_NUPB * NUP1_ABCn / VolNPC; %NUP1_ABCn –> NUP1 + ABCy
NUP1_ABCFn_Dissociation = koff_NUPB * NUP1_ABCFn / VolNPC; %NUP1_ABCFn –> NUP1 + ABCFy
NUP1_RdNn_Dissociation = koff_NUPN * NUP1_RdNn / VolNPC; %NUP1_RdNn –> NUP1 + RdNy
NUP1_RtBn_Dissociation = koff_NUPB * NUP1_RtBn / VolNPC; %NUP1_RtBn –> NUP1 + RtBy
NUP1_ABCBP80n_Dissociation = koff_NUPB * NUP1_ABCBP80n / VolNPC; %NUP1_ABCBP80n –> NUP1 + ABCBP80y
NUP1_XPO1Nn_Dissociation = koff_NUPB * NUP1_XPO1Nn / VolNPC; %NUP1_XPO1Nn –> NUP1 + XPO1Ny
NUP1_ABNNESn_Dissociation = koff_NUPB * NUP1_ABNNESn / VolNPC; %NUP1_ABNNESn –> NUP1 + ABNNESy
NUP1_XPO1Cn_Dissociation = koff_NUPB * NUP1_XPO1Cn / VolNPC; %NUP1_XPO1Cn –> NUP1 + XPO1Cy
NUP1_ABCNESn_Dissociation = koff_NUPB * NUP1_ABCNESn / VolNPC; %NUP1_ABCNESn –> NUP1 + ABCNESy
NUP1_XPO1In_Dissociation = koff_NUPB * NUP1_XPO1In / VolNPC; %NUP1_XPO1In –> NUP1 + XPO1Iy
NUP1_BINESn_Dissociation = koff_NUPB * NUP1_BINESn / VolNPC; %NUP1_BINESn –> NUP1 + BINESy
%Reactions involving RanBP3 transport; assume same properties as NLS cargo
ABy_P3y_Binding = kon1_ABC * (ABy / VolCyto) * (P3y / VolCyto) ...
    - koff1_ABC * ABP3ytemp/ VolCyto; %ABy + P3y <–> ABP3ytemp
ABP3ytemp_Conversion = kon2_ABC * ABP3ytemp / VolCyto ...
    - koff2_ABC * ABP3y / VolCyto; %ABP3ytemp <–> ABP3y
ABP3y_PeriPerm = P_ABC * ABP3y / VolCyto - P_ABC * ABP3yperi / VolNPC; %ABP3y <–> ABP3yperi
ABP3yperi_NUP1_Binding = kon_NUPB * (ABP3yperi / VolNPC) * (NUP1 / VolNPC); %ABP3yperi + NUP1 –> NUP1_ABP3y
NUP1_ABP3y_Dissociation = koff_NUPB * NUP1_ABP3y / VolNPC; %NUP1_ABP3y –> NUP1 + ABP3n
ABP3_NucPerm = P_ABC * ABP3y / VolCyto - P_ABC * ABP3n / VolNuc; %ABP3y <–> ABP3n
ABP3n_PeriPerm = P_ABC * ABP3n / VolNuc - P_ABC * ABP3nperi / VolNPC; %ABP3n <–> ABP3nperi
ABP3nperi_NUP1_Binding = kon_NUPB * (ABP3nperi / VolNPC) * (NUP1 / VolNPC); %ABP3nperi + NUP1 –> NUP1_ABP3n
NUP1_ABP3n_Dissociation = koff_NUPB * NUP1_ABP3n / VolNPC; %NUP1_ABP3n –> NUP1 + ABP3y
ABP3n_Rt_Exchange = kon_ABCRt * (ABP3n / VolNuc) * (Rt / VolNuc) ...
    - koff_ABCRt * (AP3n / VolNuc) * (RtBn / VolNuc); %ABP3n + Rt <–> AP3n + RtBn
RtSn_AP3n_Exchange = kon_RtSAC * (RtSn / VolNuc) * (AP3n / VolNuc) ...
    - koff_RtSAC * (RtSAn / VolNuc) * (P3n / VolNuc); %RtSn + AP3n <–> RtSAn + P3n
Sn_Rt_AP3n_Exchange = kon_RtSAC * (Sn / VolNuc) * (Rt / VolNuc) * (AP3n / VolNuc) * 1e-6 ...
    - koff_RtSAC * (RtSAn / VolNuc) * (P3n / VolNuc); %Sn + Rt + AP3n <–> RtSAn + P3n; can't use kon_RtSAC here but Jarnac code does anyways; rescaling rate parameter 
%Reactions involving secondary NUP that does not transport NLS-CBP80-IBB-like cargo
Byperi_NUP2_Binding = kon_NUPB * (Byperi / VolNPC) * (NUP2 / VolNPC); %Byperi + NUP2 –>  NUP2_By
Nyperi_NUP2_Binding = kon_NUPN * (Nyperi / VolNPC) * (NUP2 / VolNPC); %Nyperi + NUP2 –>  NUP2_Ny
Syperi_NUP2_Binding  = kon_NUPS * (Syperi / VolNPC) * (NUP2 / VolNPC); %Syperi + NUP2 –>  NUP2_Sy
Myperi_NUP2_Binding = kon_NUPM * (Myperi / VolNPC) * (NUP2 / VolNPC); %Myperi + NUP2 –>  NUP2_My
RtSAyperi_NUP2_Binding = kon_NUPS * (RtSAyperi / VolNPC) * (NUP2 / VolNPC); %RtSAyperi + NUP2 –>  NUP2_RtSAy
AByperi_NUP2_Binding = kon_NUPB * (AByperi / VolNPC) * (NUP2 / VolNPC); %AByperi + NUP2 –>  NUP2_ABy
RdNyperi_NUP2_Binding = kon_NUPN * (RdNyperi / VolNPC) * (NUP2 / VolNPC); %RdNyperi + NUP2 –>  NUP2_RdNy
RtByperi_NUP2_Binding = kon_NUPB * (RtByperi / VolNPC) * (NUP2 / VolNPC); %RtByperi + NUP2 –>  NUP2_RtBy
NUP2_By_Dissociation = koff_NUPB * NUP2_By / VolNPC; %NUP2_By –> NUP2 + Bn
NUP2_Ny_Dissociation = koff_NUPN * NUP2_Ny / VolNPC; %NUP2_Ny –> NUP2 + Nn
NUP2_Sy_Dissociation = koff_NUPS * NUP2_Sy / VolNPC; %NUP2_Sy –> NUP2 + Sn
NUP2_My_Dissociation = koff_NUPM * NUP2_My / VolNPC; %NUP2_My –> NUP2 + Mn
NUP2_RtSAy_Dissociation = koff_NUPS * NUP2_RtSAy / VolNPC; %NUP2_RtSAy –> NUP2 + RtSAn
NUP2_ABy_Dissociation = koff_NUPB * NUP2_ABy / VolNPC; %NUP2_ABy –> NUP2 + ABn
NUP2_RdNy_Dissociation = koff_NUPN * NUP2_RdNy / VolNPC; %NUP2_RdNy –> NUP2 + RdNn
NUP2_RtBy_Dissociation = koff_NUPB * NUP2_RtBy / VolNPC; %NUP2_RtBy –> NUP2 + RtBn
Bnperi_NUP2_Binding = kon_NUPB * (Bnperi / VolNPC) * (NUP2 / VolNPC); %Bnperi + NUP2 –>  NUP2_Bn
Nnperi_NUP2_Binding = kon_NUPN * (Nnperi / VolNPC) * (NUP2 / VolNPC); %Nnperi + NUP2 –>  NUP2_Nn
Snperi_NUP2_Binding  = kon_NUPS * (Snperi / VolNPC) * (NUP2 / VolNPC); %Snperi + NUP2 –>  NUP2_Sn
Mnperi_NUP2_Binding = kon_NUPM * (Mnperi / VolNPC) * (NUP2 / VolNPC); %Mnperi + NUP2 –>  NUP2_Mn
RtSAnperi_NUP2_Binding = kon_NUPS * (RtSAnperi / VolNPC) * (NUP2 / VolNPC); %RtSAnperi + NUP2 –>  NUP2_RtSAn
ABnperi_NUP2_Binding = kon_NUPB * (ABnperi / VolNPC) * (NUP2 / VolNPC); %ABnperi + NUP2 –>  NUP2_ABn
RdNnperi_NUP2_Binding = kon_NUPN * (RdNnperi / VolNPC) * (NUP2 / VolNPC); %RdNnperi + NUP2 –>  NUP2_RdNn
RtBnperi_NUP2_Binding = kon_NUPB * (RtBnperi / VolNPC) * (NUP2 / VolNPC); %RtBnperi + NUP2 –>  NUP2_RtBn
NUP2_Bn_Dissociation = koff_NUPB * NUP2_Bn / VolNPC; %NUP2_Bn –> NUP2 + By
NUP2_Nn_Dissociation = koff_NUPN * NUP2_Nn / VolNPC; %NUP2_Nn –> NUP2 + Ny
NUP2_Sn_Dissociation = koff_NUPS * NUP2_Sn / VolNPC; %NUP2_Sn –> NUP2 + Sy
NUP2_Mn_Dissociation = koff_NUPM * NUP2_Mn / VolNPC; %NUP2_Mn –> NUP2 + My
NUP2_RtSAn_Dissociation = koff_NUPS * NUP2_RtSAn / VolNPC; %NUP2_RtSAn –> NUP2 + RtSAy
NUP2_ABn_Dissociation = koff_NUPB * NUP2_ABn / VolNPC; %NUP2_ABn –> NUP2 + ABy
NUP2_RdNn_Dissociation = koff_NUPN * NUP2_RdNn / VolNPC; %NUP2_RdNn –> NUP2 + RdNy
NUP2_RtBn_Dissociation = koff_NUPB * NUP2_RtBn / VolNPC; %NUP2_RtBn –> NUP2 + RtBy
ABP3yperi_NUP2_Binding = kon_NUPB * (ABP3yperi / VolNPC) * (NUP2 / VolNPC); %ABP3yperi + NUP2 –> NUP2_ABP3y
NUP2_ABP3y_Dissociation = koff_NUPB * NUP2_ABP3y / VolNPC; %NUP2_ABP3y –> NUP2 + ABP3n
ABP3nperi_NUP2_Binding = kon_NUPB * (ABP3nperi / VolNPC) * (NUP2 / VolNPC); %ABP3nperi + NUP2 –> NUP2_ABP3n
NUP2_ABP3n_Dissociation = koff_NUPB * NUP2_ABP3n / VolNPC; %NUP2_ABP3n –> NUP2 + ABP3y
%Feedback for ErbB cargo
if FB
    fbwt = 0.02; %Baseline feedback weight of ErbBs on Imp-α and CAS; set to zero for no feedback simulation
    fbwt2 = 0.2; %Baseline feedback weight of Imp-α to CAS; set to zero for no feedback simulation
else
    fbwt = 0;
    fbwt2 = 0;
end
if max([CFy C2Fy CBP80y NNESy INESy CNESy]) == CFy
    ErbB_Feedback = fbwt * ((CFn / VolNuc) - (CFy / VolCyto)); %ErbBn –| An, Ay; ErbBy –> Sn, Sy
elseif max([CFy C2Fy CBP80y NNESy INESy CNESy]) == C2Fy
    ErbB_Feedback = fbwt * ((C2Fn / VolNuc) - (C2Fy / VolCyto)); %ErbBn –| An, Ay; ErbBy –> Sn, Sy
elseif max([CFy C2Fy CBP80y NNESy INESy CNESy]) == CBP80y
    ErbB_Feedback = fbwt * ((CBP80n / VolNuc) - (CBP80y / VolCyto)); %ErbBn –| An, Ay; ErbBy –> Sn, Sy
elseif max([CFy C2Fy CBP80y NNESy INESy CNESy]) == NNESy
    ErbB_Feedback = fbwt * ((NNESn / VolNuc) - (NNESy / VolCyto)); %ErbBn –| An, Ay; ErbBy –> Sn, Sy
elseif max([CFy C2Fy CBP80y NNESy INESy CNESy]) == INESy
    ErbB_Feedback = fbwt * ((INESn / VolNuc) - (INESy / VolCyto)); %ErbBn –| An, Ay; ErbBy –> Sn, Sy
elseif max([CFy C2Fy CBP80y NNESy INESy CNESy]) == CNESy
    ErbB_Feedback = fbwt * ((CNESn / VolNuc) - (CNESy / VolCyto)); %ErbBn –| An, Ay; ErbBy –> Sn, Sy
else
    ErbB_Feedback = 0;
end

%% Turn off unlimited transport of karyopherins if nuclear pores will be modeled
if NPs
    S_CytoPerm = 0;
    RtSA_CytoPerm = 0;
    N_CytoPerm = 0;
    RdN_NucPerm = 0;
    B_CytoPerm = 0;
    AB_CytoPerm = 0;
    C2B_NucPerm = 0;
    C2FB_NucPerm = 0;
    ABC_NucPerm = 0;
    ABCF_NucPerm = 0;
    ABCBP80_NucPerm = 0;
    RtB_NucPerm = 0;
    ABP3_NucPerm = 0;
end

%% Rate equations
dCdt = [- Rd_NucPerm - Rd_Ny_Binding + RtPy_Disassembly + RtPBy_Disassembly + RtNRPy_Disassembly ...
        + Rt_Hydrolysis; %Rd
    - RtSAy_P_Disassembly + RtPy_Disassembly - RtBy_P_Binding + RtPBy_Disassembly - RtNRy_P_Binding ...
        + RtNRPy_Disassembly - Rt_P_Binding - XPO1y_P_Disassembly - XPO1Ny_P_Disassembly ...
        - XPO1Cy_P_Disassembly - XPO1Iy_P_Disassembly; %P
    0; %RanGap
    - RCC1_Rdn_Binding - RCC1_RdNn_Exchange + RCC1_Rt_Unbinding; %RCC1
    N_CytoPerm - Rd_Ny_Binding - Ny_PeriPerm + NUP1_Nn_Dissociation + NUP2_Nn_Dissociation; %Ny
    A_CytoPerm - Ay_By_Binding + RtSAy_P_Disassembly - RtPBy_Ay_Exchange - Sy_Rty_Ay_Binding ...
        - ErbB_Feedback*fbwt2; %Ay
    B_CytoPerm - Ay_By_Binding + RtPBy_Disassembly - By_C2y_Binding - By_C2Fy_Binding ...
        - Rt_By_Binding - By_INESy_Binding - By_PeriPerm + NUP1_Bn_Dissociation ...
        + NUP2_Bn_Dissociation; %By
    S_CytoPerm + RtSAy_P_Disassembly - Sy_Rty_Ay_Binding - Sy_PeriPerm + NUP1_Sn_Dissociation ...
        + NUP2_Sn_Dissociation + ErbB_Feedback; %Sy
    XPO1y_P_Disassembly + XPO1Ny_P_Disassembly + XPO1Cy_P_Disassembly ...
        + XPO1Iy_P_Disassembly - My_PeriPerm + NUP1_Mn_Dissociation ...
        + NUP2_Mn_Dissociation; %My
    XPO1y_P_Disassembly + XPO1Ny_P_Disassembly + XPO1Cy_P_Disassembly + XPO1Iy_P_Disassembly ...
        - ABy_P3y_Binding; %P3y
    C_CytoPerm - ABy_Cy_Binding; %Cy
    C2_CytoPerm - By_C2y_Binding; %C2y
    - NR_NucPerm + RtNRPy_Disassembly; %NRy
    XPO1y_P_Disassembly + NES_CytoPerm; %NESy
    - Byperi_NUP1_Binding - Nyperi_NUP1_Binding - Syperi_NUP1_Binding - Myperi_NUP1_Binding ...
        - RtSAyperi_NUP1_Binding - AByperi_NUP1_Binding - C2Byperi_NUP1_Binding ...
        - C2FByperi_NUP1_Binding - ABCyperi_NUP1_Binding - ABCFyperi_NUP1_Binding ...
        - RdNyperi_NUP1_Binding - RtByperi_NUP1_Binding - ABCBP80yperi_NUP1_Binding ...
        - XPO1Nyperi_NUP1_Binding - ABNNESyperi_NUP1_Binding - XPO1Cyperi_NUP1_Binding ...
        - ABCNESyperi_NUP1_Binding - XPO1Iyperi_NUP1_Binding - BINESyperi_NUP1_Binding ...
        + NUP1_By_Dissociation + NUP1_Ny_Dissociation + NUP1_Sy_Dissociation ...
        + NUP1_My_Dissociation + NUP1_RtSAy_Dissociation + NUP1_ABy_Dissociation ...
        + NUP1_C2By_Dissociation + NUP1_C2FBy_Dissociation + NUP1_ABCy_Dissociation ...
        + NUP1_ABCFy_Dissociation + NUP1_RdNy_Dissociation + NUP1_RtBy_Dissociation ...
        + NUP1_ABCBP80y_Dissociation + NUP1_XPO1Ny_Dissociation + NUP1_ABNNESy_Dissociation ...
        + NUP1_XPO1Cy_Dissociation + NUP1_ABCNESy_Dissociation + NUP1_XPO1Iy_Dissociation ...
        + NUP1_BINESy_Dissociation - Bnperi_NUP1_Binding - Nnperi_NUP1_Binding ...
        - Snperi_NUP1_Binding - Mnperi_NUP1_Binding - RtSAnperi_NUP1_Binding ...
        - ABnperi_NUP1_Binding - C2Bnperi_NUP1_Binding - C2FBnperi_NUP1_Binding ...
        - ABCnperi_NUP1_Binding - ABCFnperi_NUP1_Binding - RdNnperi_NUP1_Binding ...
        - RtBnperi_NUP1_Binding - ABCBP80nperi_NUP1_Binding - XPO1Nnperi_NUP1_Binding ...
        - ABNNESnperi_NUP1_Binding - XPO1Cnperi_NUP1_Binding - ABCNESnperi_NUP1_Binding ...
        - XPO1Inperi_NUP1_Binding - BINESnperi_NUP1_Binding + NUP1_Bn_Dissociation ...
        + NUP1_Nn_Dissociation + NUP1_Sn_Dissociation + NUP1_Mn_Dissociation ...
        + NUP1_RtSAn_Dissociation + NUP1_ABn_Dissociation + NUP1_C2Bn_Dissociation ...
        + NUP1_C2FBn_Dissociation + NUP1_ABCn_Dissociation + NUP1_ABCFn_Dissociation ...
        + NUP1_RdNn_Dissociation + NUP1_RtBn_Dissociation + NUP1_ABCBP80n_Dissociation ...
        + NUP1_XPO1Nn_Dissociation + NUP1_ABNNESn_Dissociation + NUP1_XPO1Cn_Dissociation ...
        + NUP1_ABCNESn_Dissociation + NUP1_XPO1In_Dissociation + NUP1_BINESn_Dissociation ...
        - ABP3yperi_NUP1_Binding + NUP1_ABP3y_Dissociation - ABP3nperi_NUP1_Binding ...
        + NUP1_ABP3n_Dissociation; %NUP1
    - Byperi_NUP2_Binding - Nyperi_NUP2_Binding - Syperi_NUP2_Binding - Myperi_NUP2_Binding ...
        - RtSAyperi_NUP2_Binding - AByperi_NUP2_Binding - RdNyperi_NUP2_Binding ...
        - RtByperi_NUP2_Binding + NUP2_By_Dissociation + NUP2_Ny_Dissociation ...
        + NUP2_Sy_Dissociation + NUP2_My_Dissociation + NUP2_RtSAy_Dissociation ...
        + NUP2_ABy_Dissociation + NUP2_RdNy_Dissociation + NUP2_RtBy_Dissociation ...
        - Bnperi_NUP2_Binding - Nnperi_NUP2_Binding - Snperi_NUP2_Binding ...
        - Mnperi_NUP2_Binding - RtSAnperi_NUP2_Binding - ABnperi_NUP2_Binding ...
        - RdNnperi_NUP2_Binding - RtBnperi_NUP2_Binding + NUP2_Bn_Dissociation ...
        + NUP2_Nn_Dissociation + NUP2_Sn_Dissociation + NUP2_Mn_Dissociation ...
        + NUP2_RtSAn_Dissociation + NUP2_ABn_Dissociation + NUP2_RdNn_Dissociation ...
        + NUP2_RtBn_Dissociation - ABP3yperi_NUP2_Binding + NUP2_ABP3y_Dissociation ...
        - ABP3nperi_NUP2_Binding + NUP2_ABP3n_Dissociation; %NUP2
    - Rt_CytoPerm + RCC1_Rt_Unbinding - ABCn_Rt_Exchange - Rt_Sn_Binding - Sn_Rt_ACn_Exchange ...
        - ABCFn_Rt_Exchange - Sn_Rt_ACFn_Exchange - C2Bn_Rt_Exchange - C2FBn_Rt_Exchange ...
        - Rt_Bn_Binding - ABn_Rt_Exchange - Rt_NRn_Binding - Rt_P3n_Binding ...
        - Rt_MP3n_Binding - Rt_MP3NESn_Binding - Rt_Mn_Binding - Sn_Rt_ACBP80n_Exchange ...
        - ABCBP80n_Rt_Exchange - Sn_Rt_ANNESn_Exchange - ABNNESn_Rt_Exchange ...
        - Rt_MP3NNESn_Binding - Sn_Rt_ACNESn_Exchange - ABCNESn_Rt_Exchange ...
        - Rt_MP3CNESn_Binding - BINESn_Rt_Exchange - Rt_MP3INESn_Binding ...
        - ABP3n_Rt_Exchange - Sn_Rt_AP3n_Exchange; %Rt
    RdN_NucPerm - RCC1_RdNn_Exchange + NUP1_RdNy_Dissociation - RdNn_PeriPerm ...
        + NUP2_RdNy_Dissociation; %RdNn
    Rd_NucPerm - RCC1_Rdn_Binding; %Rdn
    - N_CytoPerm + RCC1_RdNn_Exchange + NUP1_Ny_Dissociation - Nn_PeriPerm ...
        + NUP2_Ny_Dissociation; %Nn
    - C2_CytoPerm + C2Bn_Rt_Exchange; %C2n
    - S_CytoPerm - Rt_Sn_Binding - Sn_Rt_ACn_Exchange - Sn_Rt_ACFn_Exchange ...
        - Sn_Rt_ACBP80n_Exchange - Sn_Rt_ANNESn_Exchange - Sn_Rt_ACNESn_Exchange ...
        + NUP1_Sy_Dissociation - Sn_PeriPerm + NUP2_Sy_Dissociation ...
        - Sn_Rt_AP3n_Exchange + ErbB_Feedback; %Sn
    - C_CytoPerm + RtSn_ACn_Exchange + Sn_Rt_ACn_Exchange; %Cn
    - A_CytoPerm - An_Bn_Binding + ABn_Rt_Exchange - ErbB_Feedback*fbwt2; %An
    - B_CytoPerm - An_Bn_Binding - Rt_Bn_Binding + NUP1_By_Dissociation ...
        - Bn_PeriPerm + NUP2_By_Dissociation; %Bn
    - AB_CytoPerm + ABntemp_Conversion - ABn_Rt_Exchange + NUP1_ABy_Dissociation ...
        - ABn_PeriPerm + NUP2_ABy_Dissociation; %ABn
    ABC_NucPerm - ABCn_Rt_Exchange + NUP1_ABCy_Dissociation - ABCn_PeriPerm; %ABCn
    RtB_NucPerm + ABCn_Rt_Exchange + ABCFn_Rt_Exchange + C2Bn_Rt_Exchange ...
        + C2FBn_Rt_Exchange + Rt_Bn_Binding + ABn_Rt_Exchange + ABCBP80n_Rt_Exchange ...
        + ABNNESn_Rt_Exchange + ABCNESn_Rt_Exchange + BINESn_Rt_Exchange ...
        + NUP1_RtBy_Dissociation - RtBn_PeriPerm + NUP2_RtBy_Dissociation ...
        + ABP3n_Rt_Exchange; %RtBn
    - RtSA_CytoPerm + RtSn_ACn_Exchange + Sn_Rt_ACn_Exchange + RtSn_ACFn_Exchange ...
        + Sn_Rt_ACFn_Exchange + RtSn_ACBP80n_Exchange + Sn_Rt_ACBP80n_Exchange ...
        + RtSn_ANNESn_Exchange + Sn_Rt_ANNESn_Exchange + RtSn_ACNESn_Exchange ...
        + Sn_Rt_ACNESn_Exchange + NUP1_RtSAy_Dissociation - RtSAn_PeriPerm ...
        + NUP2_RtSAy_Dissociation + RtSn_AP3n_Exchange + Sn_Rt_AP3n_Exchange; %RtSAn
    C2B_NucPerm - C2Bn_Rt_Exchange + NUP1_C2By_Dissociation - C2Bn_PeriPerm; %C2Bn
    ABCF_NucPerm - ABCFn_Rt_Exchange + NUP1_ABCFy_Dissociation - ABCFn_PeriPerm; %ABCFn
    ABCFn_Rt_Exchange - RtSn_ACFn_Exchange - Sn_Rt_ACFn_Exchange; %ACFn
    - CF_CytoPerm + RtSn_ACFn_Exchange + Sn_Rt_ACFn_Exchange; %CFn
    - C2F_CytoPerm + C2FBn_Rt_Exchange; %C2Fn
    C2FB_NucPerm - C2FBn_Rt_Exchange + NUP1_C2FBy_Dissociation - C2FBn_PeriPerm; %C2FBn
    RCC1_Rdn_Binding + RCC1_RdNn_Exchange - RCC1R_GDP_Unbinding; %RCC1Rd
    RCC1R_GDP_Unbinding - RCC1R_GTP_Binding; %RCC1R
    RCC1R_GTP_Binding - RCC1_Rt_Unbinding; %RCC1Rt
    An_Bn_Binding - ABntemp_Conversion; %ABntemp
    NR_NucPerm - Rt_NRn_Binding; %NRn
    - RtNR_CytoPerm + Rt_NRn_Binding; %RtNRn
    ABCn_Rt_Exchange - RtSn_ACn_Exchange - Sn_Rt_ACn_Exchange; %ACn
    Rt_Sn_Binding - RtSn_ACn_Exchange - RtSn_ACFn_Exchange - RtSn_ACBP80n_Exchange ...
        - RtSn_ANNESn_Exchange - RtSn_ACNESn_Exchange - RtSn_AP3n_Exchange; %RtSn
    - RdN_NucPerm + Rd_Ny_Binding - RdNy_PeriPerm + NUP1_RdNn_Dissociation ...
        + NUP2_RdNn_Dissociation; %RdNy
    AB_CytoPerm + AyBytemp_Conversion - ABy_Cy_Binding - ABy_CFy_Binding ...
        + RtPBy_Ay_Exchange - ABy_CBP80y_Binding - ABy_NNESy_Binding - ABy_CNESy_Binding ...
        - ABy_PeriPerm + NUP1_ABn_Dissociation - ABy_P3y_Binding + NUP2_ABn_Dissociation; %ABy
    - ABC_NucPerm + ABCytemp_Conversion - ABCy_PeriPerm + NUP1_ABCn_Dissociation; %ABCy
    - RtB_NucPerm - RtBy_P_Binding + Rt_By_Binding - RtBy_PeriPerm + NUP1_RtBn_Dissociation ...
        + NUP2_RtBn_Dissociation; %RtBy
    RtSA_CytoPerm - RtSAy_P_Disassembly + Sy_Rty_Ay_Binding - RtSAy_PeriPerm ...
        + NUP1_RtSAn_Dissociation + NUP2_RtSAn_Dissociation; %RtSAy
    - C2B_NucPerm + C2Bytemp_Conversion - C2By_PeriPerm + NUP1_C2Bn_Dissociation; %C2By
    RtBy_P_Binding - RtPBy_Ay_Exchange - RtPBy_Disassembly; %RtPBy
    CF_CytoPerm - ABy_CFy_Binding; %CFy
    - ABCF_NucPerm + ABCFytemp_Conversion - ABCFy_PeriPerm + NUP1_ABCFn_Dissociation; %ABCFy
    C2F_CytoPerm - By_C2Fy_Binding; %C2Fy
    - C2FB_NucPerm + C2FBytemp_Conversion - C2FBy_PeriPerm + NUP1_C2FBn_Dissociation; %C2FBy
    ABy_Cy_Binding - ABCytemp_Conversion; %ABCytemp
    ABy_CFy_Binding - ABCFytemp_Conversion; %ABCFytemp
    By_C2y_Binding - C2Bytemp_Conversion; %C2Bytemp
    By_C2Fy_Binding - C2FBytemp_Conversion; %C2FBytemp
    RtNR_CytoPerm - RtNRy_P_Binding; %RtNRy
    Rt_CytoPerm - Rt_By_Binding - Sy_Rty_Ay_Binding - Rt_Hydrolysis - Rt_P_Binding; %Rty
    RtSAy_P_Disassembly - RtPy_Disassembly + RtPBy_Ay_Exchange + Rt_P_Binding ...
        + XPO1y_P_Disassembly + XPO1Ny_P_Disassembly + XPO1Cy_P_Disassembly ...
        + XPO1Iy_P_Disassembly; %RtPy
    Ay_By_Binding - AyBytemp_Conversion; %ABytemp
    RtNRy_P_Binding - RtNRPy_Disassembly; %RtNRPy
    - ABy_CBP80y_Binding + CBP80_CytoPerm; %CBP80y
    NNES_CytoPerm - ABy_NNESy_Binding + XPO1Ny_P_Disassembly; %NNESy
    INES_CytoPerm - By_INESy_Binding + XPO1Iy_P_Disassembly; %INESy
    CNES_CytoPerm - ABy_CNESy_Binding + XPO1Cy_P_Disassembly; %CNESy
    - ABCBP80ytemp_Conversion + ABy_CBP80y_Binding; %ABCBP80ytemp
    - ABNNESytemp_Conversion + ABy_NNESy_Binding; %ABNNESytemp
    - ABCNESytemp_Conversion + ABy_CNESy_Binding; %ABCNESytemp
    - ABCBP80_NucPerm + ABCBP80ytemp_Conversion - ABCBP80y_PeriPerm + NUP1_ABCBP80n_Dissociation; %ABCBP80y
    ABCNESytemp_Conversion - ABCNESy_PeriPerm + NUP1_ABCNESn_Dissociation; %ABCNESy
    BINESy_Conversion - BINESy_PeriPerm + NUP1_BINESn_Dissociation; %BINESy
    By_INESy_Binding - BINESy_Conversion; %BINESytemp
    ABNNESytemp_Conversion - ABNNESy_PeriPerm + NUP1_ABNNESn_Dissociation; %ABNNESy
    XPO1_CytoPerm - XPO1y_P_Disassembly; %XPO1y
    - XPO1Ny_P_Disassembly - XPO1Ny_PeriPerm + NUP1_XPO1Nn_Dissociation; %XPO1Ny
    - XPO1Iy_P_Disassembly - XPO1Iy_PeriPerm + NUP1_XPO1In_Dissociation; %XPO1Iy
    - XPO1Cy_P_Disassembly - XPO1Cy_PeriPerm + NUP1_XPO1Cn_Dissociation; %XPO1Cy
    RtSn_ACBP80n_Exchange + Sn_Rt_ACBP80n_Exchange - CBP80_CytoPerm; %CBP80n
    - MP3n_NESn_Binding - RtMP3n_NESn_Binding - NES_CytoPerm; %NESn
    - NNES_CytoPerm + RtSn_ANNESn_Exchange + Sn_Rt_ANNESn_Exchange ...
        - MP3n_NNESn_Binding - RtMP3n_NNESn_Binding; %NNESn
    - INES_CytoPerm + BINESn_Rt_Exchange - MP3n_INESn_Binding - RtMP3n_INESn_Binding; %INESn
    - CNES_CytoPerm + RtSn_ACNESn_Exchange + Sn_Rt_ACNESn_Exchange ...
        - MP3n_CNESn_Binding - RtMP3n_CNESn_Binding; %CNESn
    - Mn_P3n_Binding - Rt_Mn_Binding - RtP3n_Mn_Binding + NUP1_My_Dissociation ...
        - Mn_PeriPerm + NUP2_My_Dissociation; %Mn
    - Rt_P3n_Binding - Mn_P3n_Binding - RtMn_P3n_Binding + RtSn_AP3n_Exchange ...
        + Sn_Rt_AP3n_Exchange; %P3n
    Mn_P3n_Binding - MP3n_NESn_Binding - Rt_MP3n_Binding - MP3n_NNESn_Binding ...
        - MP3n_CNESn_Binding - MP3n_INESn_Binding; %MP3n
    Rt_Mn_Binding - RtMn_P3n_Binding; %RtMn
    Rt_P3n_Binding - RtP3n_Mn_Binding; %RtP3n
    Rt_MP3n_Binding - RtMP3n_NESn_Binding + RtMn_P3n_Binding + RtP3n_Mn_Binding ...
        - RtMP3n_NNESn_Binding - RtMP3n_CNESn_Binding - RtMP3n_INESn_Binding; %RtMP3n
    MP3n_NESn_Binding - Rt_MP3NESn_Binding; %MP3NESn
    MP3n_NNESn_Binding - Rt_MP3NNESn_Binding; %MP3NNESn
    MP3n_INESn_Binding - Rt_MP3INESn_Binding; %MP3INESn
    MP3n_CNESn_Binding - Rt_MP3CNESn_Binding; %MP3CNESn
    - XPO1_CytoPerm + Rt_MP3NESn_Binding + RtMP3n_NESn_Binding; %XPO1n
    Rt_MP3CNESn_Binding + RtMP3n_CNESn_Binding + NUP1_XPO1Cy_Dissociation ...
        - XPO1Cn_PeriPerm; %XPO1Cn
    Rt_MP3INESn_Binding + RtMP3n_INESn_Binding + NUP1_XPO1Iy_Dissociation ...
        - XPO1In_PeriPerm; %XPO1In
    Rt_MP3NNESn_Binding + RtMP3n_NNESn_Binding + NUP1_XPO1Ny_Dissociation ...
        - XPO1Nn_PeriPerm; %XPO1Nn
    ABCBP80_NucPerm - RtSn_ACBP80n_Exchange - Sn_Rt_ACBP80n_Exchange + ABCBP80n_Rt_Exchange; %ACBP80n
    - RtSn_ACNESn_Exchange - Sn_Rt_ACNESn_Exchange + ABCNESn_Rt_Exchange; %ACNESn
    - BINESn_Rt_Exchange + NUP1_BINESy_Dissociation - BINESn_PeriPerm; %BINESn
    - RtSn_ANNESn_Exchange - Sn_Rt_ANNESn_Exchange + ABNNESn_Rt_Exchange; %ANNESn
    - ABCBP80n_Rt_Exchange + NUP1_ABCBP80y_Dissociation - ABCBP80n_PeriPerm; %ABCBP80n
    - ABCNESn_Rt_Exchange + NUP1_ABCNESy_Dissociation - ABCNESn_PeriPerm; %ABCNESn
    - ABNNESn_Rt_Exchange + NUP1_ABNNESy_Dissociation - ABNNESn_PeriPerm; %ABNNESn
    By_PeriPerm - Byperi_NUP1_Binding - Byperi_NUP2_Binding; %Byperi
    Ny_PeriPerm - Nyperi_NUP1_Binding - Nyperi_NUP2_Binding; %Nyperi
    Sy_PeriPerm - Syperi_NUP1_Binding - Syperi_NUP2_Binding; %Syperi
    My_PeriPerm - Myperi_NUP1_Binding - Myperi_NUP2_Binding; %Myperi
    C2By_PeriPerm - C2Byperi_NUP1_Binding; %C2Byperi
    C2FBy_PeriPerm - C2FByperi_NUP1_Binding; %C2FByperi
    RdNy_PeriPerm - RdNyperi_NUP1_Binding - RdNyperi_NUP2_Binding; %RdNyperi
    RtBy_PeriPerm - RtByperi_NUP1_Binding - RtByperi_NUP2_Binding; %RtByperi
    RtSAy_PeriPerm - RtSAyperi_NUP1_Binding - RtSAyperi_NUP2_Binding; %RtSAyperi
    ABy_PeriPerm - AByperi_NUP1_Binding - AByperi_NUP2_Binding; %AByperi
    ABCy_PeriPerm - ABCyperi_NUP1_Binding; %ABCyperi
    ABCFy_PeriPerm - ABCFyperi_NUP1_Binding; %ABCFyperi
    ABCBP80y_PeriPerm - ABCBP80yperi_NUP1_Binding; %ABCBP80yperi
    ABNNESy_PeriPerm - ABNNESyperi_NUP1_Binding; %ABNNESyperi
    ABCNESy_PeriPerm - ABCNESyperi_NUP1_Binding; %ABCNESyperi
    BINESy_PeriPerm - BINESyperi_NUP1_Binding; %BINESyperi
    XPO1Ny_PeriPerm - XPO1Nyperi_NUP1_Binding; %XPO1Nyperi
    XPO1Iy_PeriPerm - XPO1Iyperi_NUP1_Binding; %XPO1Iyperi
    XPO1Cy_PeriPerm - XPO1Cyperi_NUP1_Binding; %XPO1Cyperi
    Bn_PeriPerm - Bnperi_NUP1_Binding - Bnperi_NUP2_Binding; %Bnperi
    Nn_PeriPerm - Nnperi_NUP1_Binding - Nnperi_NUP2_Binding; %Nnperi
    Sn_PeriPerm - Snperi_NUP1_Binding - Snperi_NUP2_Binding; %Snperi
    Mn_PeriPerm - Mnperi_NUP1_Binding - Mnperi_NUP2_Binding; %Mnperi
    C2Bn_PeriPerm - C2Bnperi_NUP1_Binding; %C2Bnperi
    C2FBn_PeriPerm - C2FBnperi_NUP1_Binding; %C2FBnperi
    RdNn_PeriPerm - RdNnperi_NUP1_Binding - RdNnperi_NUP2_Binding; %RdNnperi
    RtBn_PeriPerm - RtBnperi_NUP1_Binding - RtBnperi_NUP2_Binding; %RtBnperi
    RtSAn_PeriPerm - RtSAnperi_NUP1_Binding - RtSAnperi_NUP2_Binding; %RtSAnperi
    ABn_PeriPerm - ABnperi_NUP1_Binding - ABnperi_NUP2_Binding; %ABnperi
    ABCn_PeriPerm - ABCnperi_NUP1_Binding; %ABCnperi
    ABCFn_PeriPerm - ABCFnperi_NUP1_Binding; %ABCFnperi
    ABCBP80n_PeriPerm - ABCBP80nperi_NUP1_Binding; %ABCBP80nperi
    ABNNESn_PeriPerm - ABNNESnperi_NUP1_Binding; %ABNNESnperi
    ABCNESn_PeriPerm - ABCNESnperi_NUP1_Binding; %ABCNESnperi
    BINESn_PeriPerm - BINESnperi_NUP1_Binding; %BINESnperi
    XPO1Nn_PeriPerm - XPO1Nnperi_NUP1_Binding; %XPO1Nnperi
    XPO1In_PeriPerm - XPO1Inperi_NUP1_Binding; %XPO1Inperi
    XPO1Cn_PeriPerm - XPO1Cnperi_NUP1_Binding; %XPO1Cnperi
    Bnperi_NUP1_Binding - NUP1_Bn_Dissociation; %NUP1_Bn
    Nnperi_NUP1_Binding - NUP1_Nn_Dissociation; %NUP1_Nn
    Snperi_NUP1_Binding - NUP1_Sn_Dissociation; %NUP1_Sn
    Mnperi_NUP1_Binding - NUP1_Mn_Dissociation; %NUP1_Mn
    C2Bnperi_NUP1_Binding - NUP1_C2Bn_Dissociation; %NUP1_C2Bn
    C2FBnperi_NUP1_Binding - NUP1_C2FBn_Dissociation; %NUP1_C2FBn
    RdNnperi_NUP1_Binding - NUP1_RdNn_Dissociation; %NUP1_RdNn
    RtBnperi_NUP1_Binding - NUP1_RtBn_Dissociation; %NUP1_RtBn
    RtSAnperi_NUP1_Binding - NUP1_RtSAn_Dissociation; %NUP1_RtSAn
    ABnperi_NUP1_Binding - NUP1_ABn_Dissociation; %NUP1_ABn
    ABCnperi_NUP1_Binding - NUP1_ABCn_Dissociation; %NUP1_ABCn
    ABCFnperi_NUP1_Binding - NUP1_ABCFn_Dissociation; %NUP1_ABCFn
    ABCBP80nperi_NUP1_Binding - NUP1_ABCBP80n_Dissociation; %NUP1_ABCBP80n
    ABNNESnperi_NUP1_Binding - NUP1_ABNNESn_Dissociation; %NUP1_ABNNESn
    ABCNESnperi_NUP1_Binding - NUP1_ABCNESn_Dissociation; %NUP1_ABCNESn
    BINESnperi_NUP1_Binding - NUP1_BINESn_Dissociation; %NUP1_BINESn
    XPO1Nnperi_NUP1_Binding - NUP1_XPO1Nn_Dissociation; %NUP1_XPO1Nn
    XPO1Inperi_NUP1_Binding - NUP1_XPO1In_Dissociation; %NUP1_XPO1In
    XPO1Cnperi_NUP1_Binding - NUP1_XPO1Cn_Dissociation; %NUP1_XPO1Cn
    Byperi_NUP1_Binding - NUP1_By_Dissociation; %NUP1_By
    Nyperi_NUP1_Binding - NUP1_Ny_Dissociation; %NUP1_Ny
    Syperi_NUP1_Binding - NUP1_Sy_Dissociation; %NUP1_Sy
    Myperi_NUP1_Binding - NUP1_My_Dissociation; %NUP1_My
    C2Byperi_NUP1_Binding - NUP1_C2By_Dissociation; %NUP1_C2By
    C2FByperi_NUP1_Binding - NUP1_C2FBy_Dissociation; %NUP1_C2FBy
    RdNyperi_NUP1_Binding - NUP1_RdNy_Dissociation; %NUP1_RdNy
    RtByperi_NUP1_Binding - NUP1_RtBy_Dissociation; %NUP1_RtBy
    RtSAyperi_NUP1_Binding - NUP1_RtSAy_Dissociation; %NUP1_RtSAy
    AByperi_NUP1_Binding - NUP1_ABy_Dissociation; %NUP1_ABy
    ABCyperi_NUP1_Binding - NUP1_ABCy_Dissociation; %NUP1_ABCy
    ABCFyperi_NUP1_Binding - NUP1_ABCFy_Dissociation; %NUP1_ABCFy
    ABCBP80yperi_NUP1_Binding - NUP1_ABCBP80y_Dissociation; %NUP1_ABCBP80y
    ABNNESyperi_NUP1_Binding - NUP1_ABNNESy_Dissociation; %NUP1_ABNNESy
    ABCNESyperi_NUP1_Binding - NUP1_ABCNESy_Dissociation; %NUP1_ABCNESy
    BINESyperi_NUP1_Binding - NUP1_BINESy_Dissociation; %NUP1_BINESy
    XPO1Nyperi_NUP1_Binding - NUP1_XPO1Ny_Dissociation; %NUP1_XPO1Ny
    XPO1Iyperi_NUP1_Binding - NUP1_XPO1Iy_Dissociation; %NUP1_XPO1Iy
    XPO1Cyperi_NUP1_Binding - NUP1_XPO1Cy_Dissociation; %NUP1_XPO1Cy
    - ABP3_NucPerm + ABP3ytemp_Conversion - ABP3y_PeriPerm ...
        + NUP1_ABP3n_Dissociation + NUP2_ABP3n_Dissociation; %ABP3y
    ABy_P3y_Binding - ABP3ytemp_Conversion; %ABP3ytemp
    ABP3y_PeriPerm - ABP3yperi_NUP1_Binding - ABP3yperi_NUP2_Binding; %ABP3yperi
    ABP3yperi_NUP1_Binding - NUP1_ABP3y_Dissociation; %NUP1_ABP3y
    ABP3_NucPerm + NUP1_ABP3y_Dissociation + NUP2_ABP3y_Dissociation ...
        - ABP3n_PeriPerm - ABP3n_Rt_Exchange; %ABP3n
    ABP3n_PeriPerm - ABP3nperi_NUP1_Binding - ABP3nperi_NUP2_Binding; %ABP3nperi
    ABP3nperi_NUP1_Binding - NUP1_ABP3n_Dissociation; %NUP1_ABP3n
    ABP3n_Rt_Exchange - RtSn_AP3n_Exchange - Sn_Rt_AP3n_Exchange; %AP3n
    Bnperi_NUP2_Binding - NUP2_Bn_Dissociation; %NUP2_Bn
    Nnperi_NUP2_Binding - NUP2_Nn_Dissociation; %NUP2_Nn
    Snperi_NUP2_Binding - NUP2_Sn_Dissociation; %NUP2_Sn
    Mnperi_NUP2_Binding - NUP2_Mn_Dissociation; %NUP2_Mn
    RdNnperi_NUP2_Binding - NUP2_RdNn_Dissociation; % NUP2_RdNn
    RtBnperi_NUP2_Binding - NUP2_RtBn_Dissociation; %NUP2_RtBn
    RtSAnperi_NUP2_Binding - NUP2_RtSAn_Dissociation; %NUP2_RtSAn
    ABnperi_NUP2_Binding - NUP2_ABn_Dissociation; %NUP2_ABn
    Byperi_NUP2_Binding - NUP2_By_Dissociation; % NUP2_By
    Nyperi_NUP2_Binding - NUP2_Ny_Dissociation; %NUP2_Ny
    Syperi_NUP2_Binding - NUP2_Sy_Dissociation; %NUP2_Sy
    Myperi_NUP2_Binding - NUP2_My_Dissociation; %NUP2_My
    RdNyperi_NUP2_Binding - NUP2_RdNy_Dissociation; %NUP2_RdNy
    RtByperi_NUP2_Binding - NUP2_RtBy_Dissociation; %NUP2_RtBy
    RtSAyperi_NUP2_Binding - NUP2_RtSAy_Dissociation; %NUP2_RtSAy
    AByperi_NUP2_Binding - NUP2_ABy_Dissociation; %NUP2_ABy
    ABP3yperi_NUP2_Binding - NUP2_ABP3y_Dissociation; %NUP2_ABP3y
    ABP3nperi_NUP2_Binding - NUP2_ABP3n_Dissociation]; %NUP2_ABP3n
end