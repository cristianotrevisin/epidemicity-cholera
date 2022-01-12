function [Rt, et, ratio] = diagnosis_withB(y, scenario_ocv, scenario_npi, vec1, vec2, cati_input)
    Rt = zeros(1,size(y,1));
    et = zeros(1,size(y,1));
    
    
    %Fitting Parameters
    p.theta= vec1(1);
    p.m= vec1(2);
    p.D= vec1(3);
    p.phi= vec1(4);
    p.rho = vec1(5);
    p.sigma= vec1(6);
    p.muB= vec1(7);
    p.beta0= vec1(8);
    p.psi = vec1(9);
    p.t0= round(vec1(10));
    p.r = vec1(11);

    p.b1 = vec2(1);
    p.b2 = vec2(2);
    p.t1 = round(vec2(3));
    p.t2 = round(vec2(4));
    
    %pre-defined parameters
    p.gamma=0.2;               %rate at which people recover from cholera (day^-1)
    p.mu=1/(61.4*365);         %population natality and mortality rate (day^-1)
    p.alpha=-log(0.98)*p.gamma;  %mortality rate due to cholera (day^-1)
    
    load ../data/geodata nnodes POPnodes WS_dept dist_road

    %mobility (gravity model)
    fluxes=exp(-dist_road/p.D)*diag(POPnodes);
    fluxes(sub2ind([nnodes nnodes],1:nnodes,1:nnodes))=0;
    fluxes=fluxes./repmat(sum(fluxes,2),1,nnodes);
    %to dept
    fluxes=WS_dept'*fluxes*WS_dept;
    fluxes(sub2ind([10 10],1:10,1:10))=0;
    fluxes=fluxes./repmat(sum(fluxes,2),1,10);
    POPnodes = POPnodes' * WS_dept;
    POPnodes = POPnodes';
    
    %set time
    t_initial=datenum('20.10.2010','dd.mm.yyyy');
    if scenario_npi == 2 || scenario_npi == 3 ||  scenario_npi == 4 ||  scenario_npi == 9
        t_final=datenum('01.07.2017','dd.mm.yyyy');
    else
        t_final=datenum('13.01.2019','dd.mm.yyyy');
    end
    time_model=t_initial:t_final;
    
    %extract rainfall data
    load ../data/precipitation rainfall_day date_list
    index_rainfall=find(date_list==time_model(1)):find(date_list==time_model(end));
    rainfall_day=rainfall_day(:,index_rainfall);
    
   % SCENARIOS
   if (nargin == 6) && (scenario_npi == 5 || scenario_npi == 6 ||scenario_npi == 7) 
    [ocv, ~] = generate_scenario(scenario_ocv, scenario_npi, time_model);
    cati = cati_input;
elseif (nargin ~= 6) && (scenario_npi == 5 || scenario_npi == 6 ||scenario_npi == 7)
    error("NPIs required for this scenario.")
else
    [ocv, cati] = generate_scenario(scenario_ocv, scenario_npi, time_model);
end
    
    C=y(:,6:14:end)';

    S0 = y(:,1:14:end)';        
    S1 = y(:,7:14:end)';
    S2 = y(:,11:14:end)';
    B  = y(:,5:14:end)';
    

    H = diag(POPnodes);
    
    
    % Housekeeping
    n=10;
    u=ones(1,n); 
    u1n=1:n; 
    U=sparse(u1n,u1n,u,n,n); 
    Z=sparse(zeros(n)); 
    
    for i = 1:size(y,1)
        
    
        % Loading NPI
        if i-p.t1 < 0
            keff1 = i;
        else
            keff1 = p.t1;
        end
            
        Yeff1 = ((cati(i,:) - cati(i - keff1 +1,:))'./POPnodes).^p.b1;
        
        if i-p.t2 < 0
            keff2 = i;
        else
            keff2 = p.t2;
        end
            
        Yeff2 = ((cati(i,:) - cati(i - keff2 +1,:))'./POPnodes).^p.b2;
        
        % Loading PPV
        if i<p.t0+1
            Ceff = C(:,i);
        else
            Ceff = C(:,i)-C(:,i-p.t0);
        end

        eta1 = diag(ocv.eta_1d(:,i)); eta2 = diag(ocv.eta_2d(:,i));

        % Loading EXP, SHD
        beta_t = p.beta0*exp(-Ceff./POPnodes/p.psi-Yeff1);        
        theta_t = p.theta*(1+p.phi*rainfall_day(:,i)).*exp(-Yeff2);
        
        betaM = diag(beta_t);
        thetaM = diag(theta_t);
        
        Bt = B(:,i);
        DD = diag(1./(1+Bt).^2);
        
        % Normalise populations
        s0 = diag(S0(:,i))*inv(H);
        
        % Get derived force of infections
%         TS  = ((1-p.m)*U + p.m*fluxes)*s0*betaM*DD;
TS  = s0*betaM*DD;
        % Decay of infectious
        
        PHI = p.gamma + p.alpha + p.mu;
        
        % JACOBIAN MATRIX
        
        Jac = [-PHI*U       p.sigma*TS;
                thetaM     -p.muB*U];
        
        HR = 0.5*(Jac+Jac');
        
        % TRANSMISSION MATRIX
        TR = [Z     p.sigma*TS;
              thetaM     Z];
        
        % TRANSITION MATRIX
        SIG = [-PHI*U       Z;
                Z     -p.muB*U];
        
        % DYNAMIC REPRODUCTION NUMBER
        
        NGM = -TR*inv(SIG);
        Rt(i) = max(real(eig(full(NGM))));
        
        % EPIDEMICITY
        
        et(i) = eigs(HR,1,'largestreal');

        ratio(i) = (eigs(0.5*(TR+TR'),1,'largestreal') - PHI)/(eigs(TR,1,'largestreal')/PHI);
    end
end