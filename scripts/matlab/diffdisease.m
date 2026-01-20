%% This code generates ODE simulations across different disease parameters
%% (Different reproduction numbers R0 and primary failure rates ε)
%% Breakthrough infections peak at intermediate vaccine coverage

% Note: 12/11/25:- I use the ode23t solver because the low value of m makes
% the problem stiff. I tested ode45, ode15s, ode23s and ode23t. ode45 is
% accurate but takes very long. ode15s and ode23s give blatently incorrect
% solutions. ode23t gives the correct solution as well as takes less time.


clear all; close all;

%% Parameter values
R0 = [4, 15]; %Reproduction numbers
Vareps = [0.03, 0.2]; %vaccine failure probabilities

% Read estimated assortativity values from R output (assortativity_estim.R)
phi_file = '../../output/tables/phi_estimates.csv';
if isfile(phi_file)
    phi_table = readtable(phi_file);
    phi_mean = phi_table.phi_mean;
    phi_lower = phi_table.phi_lower;
    phi_upper = phi_table.phi_upper;
    fprintf('Loaded phi estimates from %s:\n', phi_file);
    fprintf('  phi_mean:  %.3f\n', phi_mean);
    fprintf('  phi_lower: %.3f\n', phi_lower);
    fprintf('  phi_upper: %.3f\n', phi_upper);
    PHI = unique([0 phi_lower phi_mean phi_upper .9 .98]);
else
    warning('phi_estimates.csv not found. Using default PHI values.');
    warning('Run assortativity_estim.R first to generate phi estimates.');
    PHI = [0 .3 .6 .9 .98]; % Default values if file not found
end
fprintf('PHI values: %s\n', mat2str(PHI, 3));

P = [0:.01:0.9 .905:.005:1]; % Range of vaccine coverage values
Tf = 500000; % Final timepoint of the epidemic (in days)

%% Table to store steady state properties of the epidemic for different values of assortatitvity and vaccine coverage
VariableNames = {'R0','vaccine_failure','coverage','assortativity','SU','SV','IU','IV','RU','RV','fV','incidence_U','incidence_V'};
VariableTypes = {'double','double','double','double','double','double','double','double','double','double','double','double','double'};

SteadyStates = table('Size',[0,numel(VariableNames)], 'VariableTypes',VariableTypes,'VariableNames',VariableNames); 


%% Simulate Epidemics
poolobj = parpool(8);

tic

parfor ii = 1:length(P)*length(PHI)*length(R0)*length(Vareps)
    
    % Model parameters
    [i,j,k,l] = ind2sub([length(PHI), length(P), length(R0), length(Vareps)],ii);
    
    Params = model_parameters("R0",R0,"phi",PHI(i),"p",P(j),"R0",R0(k),"vareps",Vareps(l));

    % Initial conditions (mixed initial infections, in proportion to coverage)
    N_U = (1-Params.p)*Params.NT;
    N_V = Params.p*Params.NT;
    I_U = 100*(1-Params.p);
    I_V = 100*Params.p;
    R_U = 0;
    R_V = 0;

    y0 = [N_U - I_U; N_V - I_V; I_U; I_V; R_U; R_V];


    % Solve ODE system
    [t,trajectory] = ode23t(@(t,y)SIR_vaccinated_assortativity(t,y,Params),[0 Tf],y0, odeset("NonNegative",1:6,"RelTol",1e-8));

    y = trajectory(end,:);
    TotalPopulation  = sum(y);
    if abs(TotalPopulation - Params.NT) > 1e-3
        error(sprintf("Total pop not conserved, P = %.2f, Phi = %.2f",P(j), PHI(i)));
    else
        incidence_u = (Params.gamma+Params.m)*y(3)*365*1e5/TotalPopulation; % Annual disease incidence in the unvaccinated population at steady state per 100K people
        incidence_v = (Params.gamma+Params.m)*y(4)*365*1e5/TotalPopulation; % Annual disease incidence in the vaccinated population at steady state per 100K people
        SteadyStates{ii,:} = [R0(k) Params.vareps Params.p Params.phi y y(4)/(y(4)+y(3)) incidence_u incidence_v];
    end

    

end
toc

delete(poolobj);

%% Store data
% Save to the standard location for R scripts to read
output_dir = '../../data/generated';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Main output file (used by R scripts)
writetable(SteadyStates, fullfile(output_dir, 'diffdisease-ode-output.csv'));
fprintf('Saved: %s/diffdisease-ode-output.csv\n', output_dir);

% Also save timestamped backup for reproducibility
timestamp = string(datetime('now'));
timestamp = strrep(timestamp,' ','_');
timestamp = strrep(timestamp,':','-');
backup_csv = sprintf("diffdisease_%s.csv", timestamp);
writetable(SteadyStates, backup_csv);
backup_mat = sprintf("diffdisease_%s.mat", timestamp);
save(backup_mat);
fprintf('Backup saved: %s, %s\n', backup_csv, backup_mat);

%% Generate plots

for ii = 1:2
    for jj = 1:2
        
        SteadyStates0 = SteadyStates(SteadyStates.R0 == R0(ii),:);
        SteadyStates1 = SteadyStates0(SteadyStates0.vaccine_failure == Vareps(jj),:);

        f = figure;
        tiles = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");
        
        linecolors = reshape(linspace(0,0.8,length(PHI)),[],1)*[1 1 1];
        
        for i = 1:length(PHI)
        
            xdata = SteadyStates1{SteadyStates1.assortativity == PHI(i),"coverage"};
            y1data = SteadyStates1{SteadyStates1.assortativity == PHI(i),"incidence_V"};
            y2data = SteadyStates1{SteadyStates1.assortativity == PHI(i),"incidence_U"};
        
            nexttile(1);
            
            semilogy(xdata,y1data,"LineWidth",2,"Color",linecolors(end+1-i,:),"LineStyle",'-');
            hold on;
            semilogy(xdata,y2data,"LineWidth",2,"Color",linecolors(end+1-i,:),"LineStyle",'--');
            
            nexttile(2);
            plot(xdata,y1data,"LineWidth",2,"Color",linecolors(end+1-i,:),"LineStyle",'-','DisplayName',sprintf("$\\Phi = %.2f$",PHI(i)));
            hold on;
            
        end
        
        for i = 1:length(PHI)
            xdata = SteadyStates1{SteadyStates1.assortativity == PHI(i),"coverage"};
            y1data = SteadyStates1{SteadyStates1.assortativity == PHI(i),"incidence_V"};
            idx = find(y1data == max(y1data));
            idx = idx(1);
            plot(xdata(idx),y1data(idx),'o',"MarkerEdgeColor",linecolors(end+1-i,:),"MarkerFaceColor",linecolors(end+1-i,:));
        end
        
        nexttile(1);
        set(gca,'Box','off','FontSize',18,'TickLabelInterpreter','latex','LineWidth',1, ...
            'YLim',[1e-1 5e4]);
        xlabel("Vaccine coverage (p)", "Interpreter","latex");
        ylabel("Infections (log-scaled, annual per 100K)", "Interpreter","latex");
        
        nexttile(2);
        set(gca,'Box','off','FontSize',18,'TickLabelInterpreter','latex','LineWidth',1);
        legend(arrayfun(@(p) sprintf("$\\Phi = %.2f$",p),PHI,"UniformOutput",false),'Interpreter','latex','FontSize',16,'Location','best');
        xlabel("Vaccine coverage (p)", "Interpreter","latex");
        ylabel("Infections (annual per 100K)", "Interpreter","latex");
        
        title(tiles,sprintf("$R_0 = %d, \\varepsilon = %.2f$",R0(ii),Vareps(jj)),"Interpreter","latex");
    end
end

        
