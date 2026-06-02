%% This code generates ODE simulations for many assortativity values
%% Generates: manyphi-ode-output.csv 
%% Used by supplemental_figs.R to create supp-underrep.pdf

% Note: 12/11/25:- I use the ode23t solver because the low value of m makes
% the problem stiff. I tested ode45, ode15s, ode23s and ode23t. ode45 is
% accurate but takes very long. ode15s and ode23s give blatently incorrect
% solutions. ode23t gives the correct solution as well as takes less time.

clear all; close all;

output_dir = '.'; 

%% Parameter values
R0 = 15;

% Many PHI values (for exploring underreporting)
% Include phi=0 as baseline (no assortativity)
PHI = [0:.01:.99];

% Range of vaccine coverage values (hard-coded from state-level values)
P_main = [.886, .916, .948, .9, .912, .932];

Tf = 500000; % Final timepoint of the epidemic (in days)

%% Table to store steady state properties
VariableNames = {'R0','vaccine_failure','coverage','phi','SU','SV','IU','IV','RU','RV','fV','incidence_U','incidence_V'};
VariableTypes = {'double','double','double','double','double','double','double','double','double','double','double','double','double'};

%% Run main simulation (standard coverage grid)
fprintf('\n=== Running main simulation (standard coverage grid) ===\n');
SteadyStates_main = table('Size',[0,numel(VariableNames)], 'VariableTypes',VariableTypes,'VariableNames',VariableNames); 
tic

for ii = 1:length(P_main)*length(PHI)
    
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

    % Solve ODE system with tight tolerances to avoid numerical noise
    opts = odeset("NonNegative",1:6,"RelTol",1e-10,"AbsTol",1e-12);
    [t,trajectory] = ode23t(@(t,y)SIR_vaccinated_assortativity(t,y,Params),[0 Tf],y0, opts);

    y = trajectory(end,:);
    TotalPopulation  = sum(y);
    if abs(TotalPopulation - Params.NT) > 1e-3
        error(sprintf("Total pop not conserved, P = %.2f, Phi = %.4f",P_main(j), PHI(i)));
    else
        incidence_u = (Params.gamma+Params.m)*y(3)*365*1e5/TotalPopulation;
        incidence_v = (Params.gamma+Params.m)*y(4)*365*1e5/TotalPopulation;
        % Set very small incidence values to 0 to avoid numerical noise
        if incidence_u < 1e-6, incidence_u = 0; end
        if incidence_v < 1e-6, incidence_v = 0; end
        SteadyStates_main{ii,:} = [R0 Params.vareps Params.p Params.phi y y(4)/(y(4)+y(3)) incidence_u incidence_v];
    end
       
    if mod(ii,100) == 0
        fprintf('Iteration = %d\n', ii);
        fprintf('Time: %.1f seconds\n', toc);
        writetable(SteadyStates_main, fullfile(output_dir, 'manyphi-ode-output.csv'));
        
    else
        % don't print
    end
end

fprintf('Main simulation completed in %.1f seconds\n', toc);

tic


%% Store data
output_dir = '../../data/generated';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Main output file
writetable(SteadyStates_main, fullfile(output_dir, 'manyphi-ode-output.csv'));
fprintf('Saved: %s/manyphi-ode-output.csv\n', output_dir);

% Timestamped backups
timestamp = string(datetime('now'));
timestamp = strrep(timestamp,' ','_');
timestamp = strrep(timestamp,':','-');
backup_csv = sprintf("manyphi_%s.csv", timestamp);
writetable(SteadyStates_main, backup_csv);
backup_mat = sprintf("manyphi_%s.mat", timestamp);
save(backup_mat);
fprintf('Backup saved: %s, %s\n', backup_csv, backup_mat);