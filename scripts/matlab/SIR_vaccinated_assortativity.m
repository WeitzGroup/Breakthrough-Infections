function dydt = SIR_vaccinated_assortativity(t,y,params)
    
    % This function contains the differential equations used to model
    % epidemic transmission with assortativity described in ____

    % Author: Tapan Goel
    % Date: 1/14/2026

    %% Input Variables:
    % t: time (double)
    % y: state variable vector (double)
    % params: model parameters (struct)
    %         contains:
    %           phi: assortativity (double)
    %           b:   birth rate (double)
    %           p:   vaccine coverage (double)
    %           bet: transmission coefficient (double)
    %           m:   per capita death rate (double)
    %           vareps: vaccine failure probability (double)
    %           gamma: per capita recovery rate (double)
  

    %% Output Variables:
    % dydt: derivative of state variable vector (double)

    %% Code:
    
    S_U = y(1); % Unvaccinated susceptible individuals
    S_V = y(2); % Vaccinated susceptible individuals
    I_U = y(3); % Unvaccinated infectious individuals
    I_V = y(4); % Vaccinated infectious individuals
    R_U = y(5); % Unvaccinated removed individuals
    R_V = y(6); % Vaccinated removed individuals

    N_U = S_U + I_U + R_U; % Total unvaccinated individuals
    N_V = S_V + I_V + R_V; % Total vaccinated individuals
    
    
    contact_U = params.phi*(I_U/N_U) + (1-params.phi)*(I_U+I_V)/(N_U+N_V); % Infectious contacts for unvaccinated individuals
    contact_V = params.phi*(I_V/N_V) + (1-params.phi)*(I_U+I_V)/(N_U+N_V); % Infectious contacts for vaccinated individuals

    dydt = zeros(6,1);

    dydt(1) = params.b*(1-params.p) - params.bet*S_U*contact_U - params.m*S_U;

    dydt(2) = params.b*params.p*params.vareps - params.bet*S_V*contact_V - params.m*S_V;

    dydt(3) = params.bet*S_U*contact_U - (params.gamma + params.m)*I_U;

    dydt(4) = params.bet*S_V*contact_V - (params.gamma + params.m)*I_V;

    dydt(5) = params.gamma*I_U - params.m*R_U;

    dydt(6) = params.b*params.p*(1-params.vareps) + params.gamma*I_V - params.m*R_V;

    

end