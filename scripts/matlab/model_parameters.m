function params = model_parameters(vals)
    
    % This function generates the parameters for the
    % SIR_vaccinated_homophily function. Input arguments are name-value
    % pairs. Each input argument is given a default value as well.

    % Author: Tapan Goel
    % Date: 1/14/2026

    %% Input Variables:
    % vals: structure to store optional input Name-Value argument pairs.
    % These include:
    %           NT: total population size (uint 32)
    %           phi: assortativity (double)           
    %           p:   vaccine coverage (double)
    %           bet: transmission coefficient (double) (in units of per day)
    %           m:   per capita death rate (double) (in units of per day)
    %           vareps: vaccine failure probability (double)
    %           gamma: per capita recovery rate (double) (in units of per day)
    %           R0: basic reproduction number (double)
  

    %% Output Variables:
    % params: structure containing model parameters

    %% Code:
    
    arguments
        vals.NT double {mustBePositive, mustBeInteger} = 1e7; 
        vals.phi double = 0.2; 
        vals.p double {mustBeLessThanOrEqual(vals.p,1)}= .97; 
        vals.vareps double {mustBeLessThanOrEqual(vals.vareps,1)}= 0.03;
        vals.gamma double  = 0.1;
        vals.m double = 1/(365*80);
        vals.R0 double = 15;
        vals.bet double = 15*(.1+1/(365*80));                
    end

    params.NT = vals.NT;
    params.phi = vals.phi;
    params.p = vals.p;
    params.vareps = vals.vareps;
    params.gamma = vals.gamma;
    params.m = vals.m;
    params.b = params.NT*params.m;

    %% Define Mesh to Solve the Problem On
    if vals.R0 ~= 15 && vals.bet ~= 15*(.1+1/(365*80))
        error("You can either specify R0 or beta, not both");
    else
        if vals.R0 ~= 15
            params.R0 = vals.R0;
            params.bet = vals.R0*(vals.m + vals.gamma);
        else
            params.bet = vals.bet;
            params.R0 = params.bet./(params.gamma + params.m);
        end
    end
end



