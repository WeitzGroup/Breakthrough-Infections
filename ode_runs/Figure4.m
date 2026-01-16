%% This code generates Figure 4 (Breakthrough infections peak at intermediate vaccine coverage) of the manuscript


% Note: 12/11/25:- I use the ode23t solver because the low value of m makes
% the problem stiff. I tested ode45, ode15s, ode23s and ode23t. ode45 is
% accurate but takes very long. ode15s and ode23s give blatently incorrect
% solutions. ode23t gives the correct solution as well as takes less time.


clear all; close all;

%% Parameter values
R0 = 15;
PHI = [0 .3 .6 .9 .98]; % Range of assortivity values
P = [0:.01:0.9 .905:.005:1]; % Range of vaccine coverage values
Tf = 500000; % Final timepoint of the epidemic (in days)

%% Table to store steady state properties of the epidemic for different values of assortatitvity and vaccine coverage
VariableNames = {'R0','vaccine_failure','coverage','assortativity','SU','SV','IU','IV','RU','RV','fV','incidence_U','incidence_V'};
VariableTypes = {'double','double','double','double','double','double','double','double','double','double','double','double','double'};

SteadyStates = table('Size',[0,numel(VariableNames)], 'VariableTypes',VariableTypes,'VariableNames',VariableNames); 


%% Simulate Epidemics
poolobj = parpool(8);

tic

parfor ii = 1:length(P)*length(PHI)
    
    % Model parameters
    [i,j] = ind2sub([length(PHI), length(P)],ii);
    
    Params = model_parameters("R0",R0,"phi",PHI(i),"p",P(j));

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
        SteadyStates{ii,:} = [R0 Params.vareps Params.p Params.phi y y(4)/(y(4)+y(3)) incidence_u incidence_v];
    end

    

end
toc

delete(poolobj);

%% Store data (include date and time to ensure files created in the same run are used)
timestamp = string(datetime('now'));
timestamp = strrep(timestamp,' ','_');
timestamp = strrep(timestamp,':','-');
filename = sprintf("Figure4_%s.csv",timestamp);
writetable(SteadyStates,filename);
filename = sprintf("Figure4_%s.mat",timestamp);
save(filename);

%% Generate plots

f = figure;
tiles = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

linecolors = reshape(linspace(0,0.8,length(PHI)),[],1)*[1 1 1];

for i = 1:length(PHI)

    xdata = SteadyStates{SteadyStates.assortativity == PHI(i),"coverage"};
    y1data = SteadyStates{SteadyStates.assortativity == PHI(i),"incidence_V"};
    y2data = SteadyStates{SteadyStates.assortativity == PHI(i),"incidence_U"};

    nexttile(1);
    
    semilogy(xdata,y1data,"LineWidth",2,"Color",linecolors(end+1-i,:),"LineStyle",'-');
    hold on;
    semilogy(xdata,y2data,"LineWidth",2,"Color",linecolors(end+1-i,:),"LineStyle",'--');
    
    nexttile(2);
    plot(xdata,y1data,"LineWidth",2,"Color",linecolors(end+1-i,:),"LineStyle",'-','DisplayName',sprintf("$\\Phi = %.2f$",PHI(i)));
    hold on;
    
end

for i = 1:length(PHI)
    xdata = SteadyStates{SteadyStates.assortativity == PHI(i),"coverage"};
    y1data = SteadyStates{SteadyStates.assortativity == PHI(i),"incidence_V"};
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



