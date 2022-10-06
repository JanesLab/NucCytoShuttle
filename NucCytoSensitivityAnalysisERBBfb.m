%Script that performs a single-parameter sensitivity analysis for initial
%conditions for six cargo spike-ins over a 32-fold range

clear
close all hidden
warning off

B2B1ct = [0.52 1.45]; %B2B1-specific nuclear and cytoplasmic volumes; measured by microscopy 
B2B1cc = [5.18; 2.69; 0.85; 0.45; 0.81; 1.52; 2.77; 5.26; 0.50; 0.15]; %B2B1-specific concentrations; measured by quantitative immunoblotting relative to HeLa 
B2B1lr = [10; 5; 6; 1]; %B2B1-specific lumped receptors (default: [10; 1; 3.6; 1])
B2B1nups = [2600; 0]; %B2B1-specific NUPs
Cargo = 1; %µM, corresponding to peak behavior among the cargo

NLSparamtable = readtable(strcat(pwd,'/NLS_variant_parameters.xlsx'));
NLSparams = [NLSparamtable{2:6,4} NLSparamtable{2:6,6} ...
      NLSparamtable{2:6,5} NLSparamtable{2:6,7}];

Scaling = 2.^linspace(-3,3,7)';

%% Sweep through initial conditions

%Sweep through cell concentrations
for i=1:height(B2B1cc)
    for j=1:height(Scaling)
        for k=1:height(NLSparams)
            B2B1cc_temp = B2B1cc;
            B2B1cc_temp(i) = B2B1cc_temp(i)*Scaling(j);
            [ot{i,j,k},cvt] = NucCytoShuttleERBBfbModel(Cargo,'CellType',B2B1ct, ...
                'CellConcs',B2B1cc_temp,'LumpedReceptors',B2B1lr,'NUPs', ...
                B2B1nups,'NLSparams',NLSparams(k,:));
        j
        end
    end
end

%% Preallocate matrices for nuclear-cytoplasmic masses
CBP80_NC = zeros(height(ot),width(ot),2);

%% Extract nuclear-cytoplasmic masses for each cargo
for i=1:height(ot)
    for j=1:width(ot)
        for k=1:height(NLSparams)
            CBP80_NC(i,j,k,1) = sum(ot{i,j,k}.CBP80{height(ot{i,j,k}.CBP80),contains(ot{i,j,k}.CBP80.Properties.VariableNames, ...
                'CBP80n')});
            CBP80_NC(i,j,k,2) = sum(ot{i,j,k}.CBP80{height(ot{i,j,k}.CBP80),contains(ot{i,j,k}.CBP80.Properties.VariableNames, ...
                'CBP80y')});
        end
    end
end

%% Plot ratio as a function of cargo concentration
xylims = [2^-4 2^4 0.3 30];
ModelSpecies = {'Rd'; 'P'; 'RanGap'; 'RCC1'; 'N'; 'A'; 'B'; ...
    'S'; 'M'; 'P3'; 'NUP1'; 'NUP2'};
for i=1:height(ot)
    subplot(3,4,i)
    loglog(Scaling,CBP80_NC(i,:,1,1)./CBP80_NC(i,:,1,2)*cvt.VolCyto/cvt.VolNuc,'g.-', ...
        'LineWidth',1,'MarkerSize',15);
    axis(xylims)
    title(ModelSpecies{i})
    xlabel('Relative abundance')
    ylabel('Nuc/cyto ratio')
    hold on
    loglog(Scaling,CBP80_NC(i,:,2,1)./CBP80_NC(i,:,2,2)*cvt.VolCyto/cvt.VolNuc,'c.-', ...
        'LineWidth',1,'MarkerSize',15);
    loglog(Scaling,CBP80_NC(i,:,3,1)./CBP80_NC(i,:,3,2)*cvt.VolCyto/cvt.VolNuc,'m.-', ...
        'LineWidth',1,'MarkerSize',15);
    loglog(Scaling,CBP80_NC(i,:,4,1)./CBP80_NC(i,:,4,2)*cvt.VolCyto/cvt.VolNuc,'y.-', ...
        'LineWidth',1,'MarkerSize',15);
    loglog(Scaling,CBP80_NC(i,:,5,1)./CBP80_NC(i,:,5,2)*cvt.VolCyto/cvt.VolNuc,'k.-', ...
        'LineWidth',1,'MarkerSize',15);
end
legend({'Kd = 100 nM' 'Kd = 30 nM' 'Kd = 10 nM' 'Kd = 3 nM' 'Kd = 1 nM'})