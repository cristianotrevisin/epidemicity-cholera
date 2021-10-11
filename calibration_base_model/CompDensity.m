function [p,log_p] = CompDensity(x,MCMCPar,Measurement,ModelName,Extra,option)
% This function computes the density of each x value

p = []; log_p = [];

% Loop over the individual parameter combinations of x
for ii = 1:size(x,1), 
    % Call model to generate simulated data
%     evalstr = ['ModPred = ',ModelName,'(x(ii,:),Extra);']; eval(evalstr);
    try
        evalstr = ['[ModPred,times_out] = ',ModelName,'(x(ii,:),Extra);']; eval(evalstr);
    catch
        ModPred = Inf;
    end
    if option == 1, % Model directly computes posterior density
        p(ii,1:2) = [ModPred ii]; log_p(ii,1) = log(p(ii,1));
    end;

    if option == 2, % Model computes output simulation
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Compute the number of measurement data
        N = size(Measurement.MeasData,1);
        % Derive the log likelihood
        log_p(ii,1) = sum(log(MCMCPar.Wb./Measurement.Sigma(:))) - MCMCPar.Cb.*(sum((abs(Err./Measurement.Sigma(:))).^(2/(1+MCMCPar.Gamma))));
        %Measurement.Sigma=100.0;
        %log_p(ii,1) = N.*log(MCMCPar.Wb./Measurement.Sigma) - MCMCPar.Cb.*(sum((abs(Err./Measurement.Sigma)).^(2/(1+MCMCPar.Gamma))))       
        % And retain in memory
        p(ii,1:2) = [log_p(ii,1) ii]; 
    end;

    if option == 3, % Model computes output simulation
        try
            Err = (Measurement.MeasData(:)-ModPred(:));
            % Derive the sum of squared error
            SSR = sum(abs(Err).^(2/(1+MCMCPar.Gamma)));         
        catch
            SSR = Inf;
        end
        % And retain in memory
        p(ii,1:2) = [-SSR ii]; log_p(ii,1) = -0.5 * SSR;
    end;

    if option == 4, % Model directly computes log posterior density
        p(ii,1:2) = [ModPred ii]; log_p(ii,1) = p(ii,1);
    end
    if option == 5, % Similar as 3, but now weights with the Measurement Sigma
        % Defime the error
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Derive the sum of squared error
        SSR = sum(abs(Err).^(2/(1+MCMCPar.Gamma)));
        % And retain in memory
        p(ii,1:2) = [-SSR ii]; log_p(ii,1) = -0.5 * SSR;
    end
    if option == 6 , % modified by DP
        if sum(isnan(ModPred(:)))==0
            ind1=find(floor(times_out)==floor(Measurement.time_data1),1,'first');
            ind2=find(floor(times_out)==floor(Measurement.time_data2),1,'last');
            
            ModPred=ModPred(:,ind1:ind2);
            % log likelihood, with a negative binomial model of observation
            %size(Measurement.MeasData)
            %size(ModPred)
            Err = (Measurement.MeasData(:)-ModPred(:));
            % Derive the sum of squared error
            SSR = sum(abs(Err).^(2));
            % And retain in memory
            p(ii,1:2) = [-SSR ii]; log_p(ii,1) = -0.5 * SSR;
        else
            p(ii,1:2) = [-Inf ii]; log_p(ii,1) = Inf;
            disp('error in convergence')  
        end

    end;
    
    if option == 7
        ind1=find(floor(times_out)==floor(Measurement.time_data1),1,'first');
        ind2=find(floor(times_out)==floor(Measurement.time_data2),1,'last');
        ModPred=ModPred(:,ind1:ind2);
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Derive the log density
        if isnan(ModPred) == 0, %&& size(Sigma,1) == 1,  % --> homoscedastic error (CK: commented out Sigma criteria since Sigma defined as 10 for now)
            log_p(ii,1) = - ( Measurement.N / 2) * log(2 * pi) - Measurement.N * log( Measurement.Sigma ) - 1/2 * Measurement.Sigma^(-2) * sum ( Err.^2 );
            
        else %isnan(ModPred) ~= 0,
            log_p(ii,1) = -1.0e+20;    %CK: probability is almost zero
        %else% --> heteroscedastic error
            %log_p(ii,1) = - ( Measurement.N / 2) * log(2 * pi) - sum ( log( Measurement.Sigma ) ) - 1/2 * sum ( ( Err./Measurement.Sigma ).^2);
        end;
        % And retain in memory
        p(ii,1) = log_p(ii,1);
        log_p(ii,2) = ii;
        p(ii,2) = ii;
    end;

end;
