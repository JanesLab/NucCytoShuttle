%Script that performs a single-parameter sensitivity analysis for initial
%conditions for six cargo spike-ins over a 32-fold range

clear
close all hidden

B2B1ct = [0.52 1.45]; %B2B1-specific nuclear and cytoplasmic volumes; measured by microscopy 
B2B1cc = [5.18; 2.69; 0.85; 0.45; 0.81; 1.52; 2.77; 5.26; 0.50; 0.15]; %B2B1-specific concentrations; measured by quantitative immunoblotting relative to HeLa 
B2B1lr = [10; 5; 6; 1]; %B2B1-specific lumped receptors (default: [10; 1; 3.6; 1])
B2B1nups = [2600; 0]; %B2B1-specific NUPs
Cargo = 1; %µM, corresponding to peak behavior among the cargo

Scaling = 2.^linspace(-3,3,7)';

%% Sweep through initial conditions

%Sweep through cell concentrations
for i=1:height(B2B1cc)
    for j=1:height(Scaling)
        B2B1cc_temp = B2B1cc;
        B2B1cc_temp(i) = B2B1cc_temp(i)*Scaling(j);
        [ot{i,j},cvt] = NucCytoShuttleModel(Cargo,'CellType',B2B1ct, ...
            'CellConcs',B2B1cc_temp,'LumpedReceptors',B2B1lr,'NUPs', ...
            B2B1nups);
        sprintf('%2d percent complete.',round(((i-1)*j+j)/(height(B2B1cc)*height(Scaling))*100))
    end
end

%Sweep through NUPs
for j=1:height(Scaling)
    B2B1nups_temp = B2B1nups;
    B2B1nups_temp = B2B1nups_temp*Scaling(j);
    [ot{height(B2B1cc)+1,j},cvt] = NucCytoShuttleModel(Cargo,'CellType',B2B1ct, ...
        'CellConcs',B2B1cc,'LumpedReceptors',B2B1lr,'NUPs', ...
        B2B1nups_temp);
end

%% Preallocate matrices for nuclear-cytoplasmic masses
NLS_NC = zeros(height(ot),width(ot),2);
IBB_NC = zeros(height(ot),width(ot),2);
CBP80_NC = zeros(height(ot),width(ot),2);
NNES_NC = zeros(height(ot),width(ot),2);
INES_NC = zeros(height(ot),width(ot),2);
CNES_NC = zeros(height(ot),width(ot),2);

%% Extract nuclear-cytoplasmic masses for each cargo
for i=1:height(ot)
    for j=1:width(ot)
        NLS_NC(i,j,1) = sum(ot{i,j}.NLS{height(ot{i,j}.NLS),contains(ot{i,j}.NLS.Properties.VariableNames, ...
            'CFn')});
        NLS_NC(i,j,2) = sum(ot{i,j}.NLS{height(ot{i,j}.NLS),contains(ot{i,j}.NLS.Properties.VariableNames, ...
            'CFy')});
        IBB_NC(i,j,1) = sum(ot{i,j}.IBB{height(ot{i,j}.IBB),contains(ot{i,j}.IBB.Properties.VariableNames, ...
            'C2Fn') | contains(ot{i,j}.IBB.Properties.VariableNames,'C2FBn')});
        IBB_NC(i,j,2) = sum(ot{i,j}.IBB{height(ot{i,j}.IBB),contains(ot{i,j}.IBB.Properties.VariableNames, ...
            'C2Fy') | contains(ot{i,j}.IBB.Properties.VariableNames,'C2FBy')});
        CBP80_NC(i,j,1) = sum(ot{i,j}.CBP80{height(ot{i,j}.CBP80),contains(ot{i,j}.CBP80.Properties.VariableNames, ...
            'CBP80n')});
        CBP80_NC(i,j,2) = sum(ot{i,j}.CBP80{height(ot{i,j}.CBP80),contains(ot{i,j}.CBP80.Properties.VariableNames, ...
            'CBP80y')});
        NNES_NC(i,j,1) = sum(ot{i,j}.NNES{height(ot{i,j}.NNES),contains(ot{i,j}.NNES.Properties.VariableNames, ...
            'NNESn') | contains(ot{i,j}.NNES.Properties.VariableNames,'XPO1Nn')});
        NNES_NC(i,j,2) = sum(ot{i,j}.NNES{height(ot{i,j}.NNES),contains(ot{i,j}.NNES.Properties.VariableNames, ...
            'NNESy') | contains(ot{i,j}.NNES.Properties.VariableNames,'XPO1Ny')});
        INES_NC(i,j,1) = sum(ot{i,j}.INES{height(ot{i,j}.INES),contains(ot{i,j}.INES.Properties.VariableNames, ...
            'INESn') | contains(ot{i,j}.NNES.Properties.VariableNames,'XPO1In')});
        INES_NC(i,j,2) = sum(ot{i,j}.INES{height(ot{i,j}.INES),contains(ot{i,j}.INES.Properties.VariableNames, ...
            'INESy') | contains(ot{i,j}.NNES.Properties.VariableNames,'XPO1Iy')});
        CNES_NC(i,j,1) = sum(ot{i,j}.CNES{height(ot{i,j}.CNES),contains(ot{i,j}.CNES.Properties.VariableNames, ...
            'CNESn') | contains(ot{i,j}.NNES.Properties.VariableNames,'XPO1Cn')});
        CNES_NC(i,j,2) = sum(ot{i,j}.CNES{height(ot{i,j}.CNES),contains(ot{i,j}.CNES.Properties.VariableNames, ...
            'CNESy') | contains(ot{i,j}.NNES.Properties.VariableNames,'XPO1Cy')});
    end
end

%% Plot ratio as a function of cargo concentration
xylims = [2^-4 2^4 0.3 30];
ModelSpecies = {'Rd'; 'P'; 'RanGap'; 'RCC1'; 'N'; 'A'; 'B'; ...
    'S'; 'M'; 'P3'; 'NUP1'; 'NUP2'};
for i=1:height(ot)
    subplot(3,4,i)
    loglog(Scaling,NLS_NC(i,:,1)./NLS_NC(i,:,2)*cvt.VolCyto/cvt.VolNuc,'c.-', ...
        'LineWidth',1,'MarkerSize',15);
    axis(xylims)
    title(ModelSpecies{i})
    xlabel('Relative abundance')
    ylabel('Nuc/cyto ratio')
    hold on
    loglog(Scaling,IBB_NC(i,:,1)./IBB_NC(i,:,2)*cvt.VolCyto/cvt.VolNuc,'m.-', ...
        'LineWidth',1,'MarkerSize',15);
    loglog(Scaling,CBP80_NC(i,:,1)./CBP80_NC(i,:,2)*cvt.VolCyto/cvt.VolNuc,'k.-', ...
        'LineWidth',1,'MarkerSize',15);
    loglog(Scaling,NNES_NC(i,:,1)./NNES_NC(i,:,2)*cvt.VolCyto/cvt.VolNuc,'c.:', ...
        'LineWidth',1,'MarkerSize',15);
    loglog(Scaling,INES_NC(i,:,1)./INES_NC(i,:,2)*cvt.VolCyto/cvt.VolNuc,'m.:', ...
        'LineWidth',1,'MarkerSize',15);
    loglog(Scaling,CNES_NC(i,:,1)./CNES_NC(i,:,2)*cvt.VolCyto/cvt.VolNuc,'k.:', ...
        'LineWidth',1,'MarkerSize',15);
end
legend({'NLS'; 'IBB'; 'CBP80'; 'NNES'; 'INES'; 'CNES'})