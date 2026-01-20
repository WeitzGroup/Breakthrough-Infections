%% This code generates ODE simulations for LOW assortativity values (near HIT)
%% Generates: lowphi-ode-output.csv and lowphi-inset-ode-output.csv
%% Used by supplemental_figs.R to create supp-lowphiHIT.pdf

% Note: 12/11/25:- I use the ode23t solver because the low value of m makes
% the problem stiff. I tested ode45, ode15s, ode23s and ode23t. ode45 is
% accurate but takes very long. ode15s and ode23s give blatently incorrect
% solutions. ode23t gives the correct solution as well as takes less time.

clear all; close all;

%% Parameter values
R0 = 15;

% Low PHI values (for exploring behavior near HIT)
PHI = [1e-04, 0.001, 0.01, 0.1];
fprintf('Low PHI values: %s\n', mat2str(PHI, 4));

% Range of vaccine coverage values (standard grid)
P_main = [0:.01:0.9 .905:.005:1];
% Finer grid near elimination threshold for inset plot
P_inset = [0.93:.001:1];

Tf = 500000; % Final timepoint of the epidemic (in days)

%% Table to store steady state properties
VariableNames = {'R0','vaccine_failure','coverage','phi','SU','SV','IU','IV','RU','RV','fV','incidence_U','incidence_V'};
VariableTypes = {'double','double','double','double','double','double','double','double','double','double','double','double','double'};

%% Run main simulation (standard coverage grid)
fprintf('\n=== Running main simulation (standard coverage grid) ===\n');
SteadyStates_main = table('Size',[0,numel(VariableNames)], 'VariableTypes',VariableTypes,'VariableNames',VariableNames); 

poolobj = parpool(8);
tic

parfor ii = 1:length(P_main)*length(PHI)
    
    % Model parameters
    [i,j] = ind2sub([length(PHI), length(P_main)],ii);
    
    Params = model_parameters("R0",R0,"phi",PHI(i),"p",P_main(j));

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
        error(sprintf("Total pop not conserved, P = %.2f, Phi = %.4f",P_main(j), PHI(i)));
    else
        incidence_u = (Params.gamma+Params.m)*y(3)*365*1e5/TotalPopulation;
        incidence_v = (Params.gamma+Params.m)*y(4)*365*1e5/TotalPopulation;
        SteadyStates_main{ii,:} = [R0 Params.vareps Params.p Params.phi y y(4)/(y(4)+y(3)) incidence_u incidence_v];
    end
end

fprintf('Main simulation completed in %.1f seconds\n', toc);

%% Run inset simulation (fine coverage grid near threshold)
fprintf('\n=== Running inset simulation (fine coverage grid) ===\n');
SteadyStates_inset = table('Size',[0,numel(VariableNames)], 'VariableTypes',VariableTypes,'VariableNames',VariableNames); 

tic

parfor ii = 1:length(P_inset)*length(PHI)
    
    % Model parameters
    [i,j] = ind2sub([length(PHI), length(P_inset)],ii);
    
    Params = model_parameters("R0",R0,"phi",PHI(i),"p",P_inset(j));

    % Initial conditions
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
        error(sprintf("Total pop not conserved, P = %.2f, Phi = %.4f",P_inset(j), PHI(i)));
    else
        incidence_u = (Params.gamma+Params.m)*y(3)*365*1e5/TotalPopulation;
        incidence_v = (Params.gamma+Params.m)*y(4)*365*1e5/TotalPopulation;
        SteadyStates_inset{ii,:} = [R0 Params.vareps Params.p Params.phi y y(4)/(y(4)+y(3)) incidence_u incidence_v];
    end
end

fprintf('Inset simulation completed in %.1f seconds\n', toc);

delete(poolobj);

%% Store data
output_dir = '../../data/generated';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Main output file
writetable(SteadyStates_main, fullfile(output_dir, 'lowphi-ode-output.csv'));
fprintf('Saved: %s/lowphi-ode-output.csv\n', output_dir);

% Inset output file (fine grid)
writetable(SteadyStates_inset, fullfile(output_dir, 'lowphi-inset-ode-output.csv'));
fprintf('Saved: %s/lowphi-inset-ode-output.csv\n', output_dir);

% Timestamped backups
timestamp = string(datetime('now'));
timestamp = strrep(timestamp,' ','_');
timestamp = strrep(timestamp,':','-');
backup_csv = sprintf("lowphi_%s.csv", timestamp);
writetable(SteadyStates_main, backup_csv);
backup_mat = sprintf("lowphi_%s.mat", timestamp);
save(backup_mat);
fprintf('Backup saved: %s, %s\n', backup_csv, backup_mat);

%% Generate verification plot
f = figure;
tiles = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

linecolors = reshape(linspace(0,0.6,length(PHI)),[],1)*[1 1 1];

for i = 1:length(PHI)
    xdata = SteadyStates_main{SteadyStates_main.phi == PHI(i),"coverage"};
    y1data = SteadyStates_main{SteadyStates_main.phi == PHI(i),"incidence_V"};
    
    nexttile(1);
    semilogy(xdata,y1data,"LineWidth",2,"Color",linecolors(end+1-i,:),"LineStyle",'-');
    hold on;
    
    nexttile(2);
    plot(xdata,y1data,"LineWidth",2,"Color",linecolors(end+1-i,:),"LineStyle",'-','DisplayName',sprintf("$\\phi = %.4f$",PHI(i)));
    hold on;
end

nexttile(1);
set(gca,'Box','off','FontSize',14,'TickLabelInterpreter','latex','LineWidth',1);
xlabel("Vaccine coverage (p)", "Interpreter","latex");
ylabel("Breakthrough infections (log, annual/100K)", "Interpreter","latex");
title("Log scale", "Interpreter","latex");

nexttile(2);
set(gca,'Box','off','FontSize',14,'TickLabelInterpreter','latex','LineWidth',1);
legend(arrayfun(@(p) sprintf("$\\phi = %.4f$",p),PHI,"UniformOutput",false),'Interpreter','latex','FontSize',12,'Location','best');
xlabel("Vaccine coverage (p)", "Interpreter","latex");
ylabel("Breakthrough infections (annual/100K)", "Interpreter","latex");
title("Linear scale", "Interpreter","latex");

sgtitle("Low $\phi$ simulations (near HIT behavior)", "Interpreter","latex","FontSize",16);

fprintf('\n=== Low phi simulation complete ===\n');
fprintf('Files generated:\n');
fprintf('  - %s/lowphi-ode-output.csv\n', output_dir);
fprintf('  - %s/lowphi-inset-ode-output.csv\n', output_dir);
