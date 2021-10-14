function [ModPred,time_model_out, y, cati, ocv, rainfall_day] = m4r(vec1, vec2)
% MODEL FOR SCENARIO ANALYSIS

    load ../data/precipitation rainfall_day date_list
    load ../data/geodata nnodes POPnodes WS_dept dist_road
    load '../data/ocv.mat' datespace rv_1d rv_2d eta_1d eta_2d
    load ../data/wash date_list_cati cati_day

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
    
    %set time
    t0=0;                    %delay between the first day of data available and the day of the onset of the epidemic (initial condition)
    t_initial=datenum('20.10.2010','dd.mm.yyyy');
    t_final=datenum('01.07.2017','dd.mm.yyyy');
    time_data=t_initial+(7-weekday(t_initial))+(0:350-1)*7; 
    
    time_model=t_initial-t0:t_final;

    %extract rainfall data
    index_rainfall=find(date_list==time_model(1)):find(date_list==time_model(end));
    rainfall_day=rainfall_day(:,index_rainfall);
   
    %mobility (gravity model)
    fluxes=exp(-dist_road/p.D)*diag(POPnodes);
    fluxes(sub2ind([nnodes nnodes],1:nnodes,1:nnodes))=0;
    fluxes=fluxes./repmat(sum(fluxes,2),1,nnodes);
    %to dept
    fluxes=WS_dept'*fluxes*WS_dept;
    fluxes(sub2ind([10 10],1:10,1:10))=0;
    fluxes=fluxes./repmat(sum(fluxes,2),1,10);
    
    %extract ocv data
    index_ocv=find(datespace==time_model(1)):find(datespace==time_model(end));
    ocv.rv_1d=rv_1d(index_ocv,:)';
    ocv.rv_2d=rv_2d(index_ocv,:)';
    ocv.eta_1d=eta_1d(index_ocv,:)';
    ocv.eta_2d=eta_2d(index_ocv,:)';
    
    %extract cati data
    index_cati=find(date_list_cati==time_model(1)):find(date_list_cati==time_model(end));
    cati_sum = cumsum(cati_day,1);
    cati = double(cati_sum(index_cati,:));
    
    % initial condition
    nnodes = 10;
    POPnodes = (POPnodes' * WS_dept)';
    
    initial_points=[1 2];
    I_initial=[1100 1000]/p.sigma;
    
    y0=zeros(14,nnodes);
    y0(1,:)=POPnodes; %whole population
    y0(1,initial_points) = y0(1,initial_points)-I_initial; %susceptible = pop minus first infected
    y0(2,initial_points) = p.sigma*I_initial; %first symptomatic infections
    y0(3,initial_points) = (1-p.sigma)*I_initial; %first removed
    y0(6,initial_points) = p.sigma*I_initial; %reported cases @ beginning
    y0(5,initial_points) = (y0(2,initial_points)+p.r*y0(3,initial_points))*p.theta/p.muB./POPnodes(initial_points)';

    %run model
    y=SIB(p,nnodes,POPnodes,fluxes,rainfall_day,cati,ocv,time_model,y0);
    cumcases_AD1_t=y(:,6:14:end)';

    %upscale from daily to weekly time scale
    index_time=zeros(length(time_data),1);
    for iii=1:length(time_data)
        index_time(iii)=find(time_model==time_data(iii));
    end

    cumcases_AD1_week=cumcases_AD1_t(:,index_time);
    cases_AD1_week=diff([zeros(size(cumcases_AD1_week,1),1) cumcases_AD1_week],1,2); %first column is 0, makes the difference between each element along the row.
    
    ModPred = cases_AD1_week(:,1:350);
    time_model_out = time_data(1:350);
    
    
    %%% ODE PART
    
    function ytop=SIB(p,nnodes,H,fluxes,rainfall_day,cati,ocv,tspan,y0)

        step_size = 7;
        for iii = 1:step_size:length(tspan)
            try
                tspan_sub = tspan(iii:iii+step_size);
            catch
                clear tspan_sub
                tspan_sub = tspan(iii:end);
            end
        
            opt=odeset('RelTol', 1e-2, 'AbsTol', 1e-3);
            [~,y]=ode45(@eqs,tspan_sub,y0,opt);
            
            for i = 1:10
                if y(end,2+14*(i-1)) < 1
                    y(end,2+14*(i-1))=0;
                end
                if y(end,3+14*(i-1)) < 1
                    y(end,3+14*(i-1))=0;
                end
                if y(end,8+14*(i-1)) < 1
                    y(end,8+14*(i-1))=0;
                end
                if y(end,9+14*(i-1)) < 1
                    y(end,9+14*(i-1))=0;
                end
                if y(end,12+14*(i-1)) < 1
                    y(end,12+14*(i-1))=0;
                end
                if y(end,13+14*(i-1)) < 1
                    y(end,13+14*(i-1))=0;
                end
            end

            try
                ytop(iii:iii+step_size,:) = y;
            catch
                ytop(iii:length(tspan),:) = y;
            end
            y0 = reshape(y(end,:),14,nnodes);
            
            clear y
        end

        function dy=eqs(t,y)
            index_t=floor(t-tspan(1))+1;

            dy=zeros(14*nnodes,1);

            temp=fluxes*(y(5:14:end)./(y(5:14:end)+1));      
            
            if index_t-p.t1 < 0
                keff1 = index_t;
            else
                keff1 = p.t1;
            end
            
            Yeff1 = ((cati(index_t,:) - cati(index_t - keff1 +1,:))'./H).^p.b1;
            
            if index_t-p.t2 < 0
                keff2 = index_t;
            else
                keff2 = p.t2;
            end
            
            Yeff2 = ((cati(index_t,:) - cati(index_t - keff2 +1,:))'./H).^p.b2;
            
            
            if index_t<p.t0+1
                Ceff = y(6:14:end);
            else
                Ceff = y(6:14:end)-ytop(index_t-p.t0,6:14:end)';
            end

            % EXPONENTIAL
            beta_t = p.beta0*exp(-Ceff./H/p.psi - Yeff1);
            
            theta_t = p.theta*(1+p.phi*rainfall_day(:,index_t)).*exp(-Yeff2);
            
            % MONO
%             beta_t = p.beta0*exp(-Ceff./H/p.psi)./(ones(10,1)+Yeff1);
%             
%             theta_t = p.theta*(1+p.phi*rainfall_day(:,index_t))./(ones(10,1)+Yeff2);

            rv1 = ocv.rv_1d(:,index_t) ./ (y(1:14:end)+y(3:14:end)+y(4:14:end));
            rv2 = zeros(10,1);
            for i = 1:10
                if y(7+14*(i-1))+y(9+14*(i-1))+y(10+14*(i-1)) < 1
                    rv2(i) = 0;
                else
                    rv2(i) = ocv.rv_2d(i,index_t) ./ (y(7+14*(i-1))+y(9+14*(i-1))+y(10+14*(i-1)));
                end
            end
            
                        
            FI=((1-p.m)*y(5:14:end)./(1+y(5:14:end))+temp*p.m).*(beta_t);
            GI = (y(2:14:end) + y(8:14:end) + y(12:14:end)) + p.r*((1-p.m)*(y(3:14:end)+y(9:14:end)+y(13:14:end)) + p.m*fluxes'*(y(3:14:end)+y(9:14:end)+y(13:14:end)));
                        
            dy(1:14:end)=p.mu*(H-y(1:14:end))- (FI+rv1).*y(1:14:end) +p.rho*y(4:14:end);        %S
            dy(2:14:end)=p.sigma.*FI.*y(1:14:end) - (p.gamma+p.alpha+p.mu)*y(2:14:end);       %I
            dy(3:14:end)=(1-p.sigma).*FI.*y(1:14:end) - (p.gamma+p.mu+rv1).*y(3:14:end);       %A
            dy(4:14:end)=p.gamma*(y(2:14:end)+y(3:14:end)) - (p.mu+p.rho+rv1).*y(4:14:end);       %R
            dy(5:14:end)=-(p.muB).*y(5:14:end)+theta_t.*GI./H;  %B
            dy(6:14:end)=p.sigma.*FI.*(y(1:14:end)+(1-ocv.eta_1d(:,index_t)).*y(7:14:end)+(1-ocv.eta_2d(:,index_t)).*y(11:14:end));   %Cum(I)  
            

            dy(7:14:end)  = rv1.*y(1:14:end)-(FI.*(1-ocv.eta_1d(:,index_t))+p.mu+rv2).*y(7:14:end)+p.rho*y(10:14:end); % V^S_1d
            dy(8:14:end)  = p.sigma.*FI.*(1-ocv.eta_1d(:,index_t)).*y(7:14:end)-(p.gamma+p.alpha+p.mu)*y(8:14:end); % V^I_1d
            dy(9:14:end)  = (1-p.sigma).*FI.*(1-ocv.eta_1d(:,index_t)).*y(7:14:end) + rv1.*y(3:14:end) - (p.gamma+p.mu+rv2).*y(9:14:end); % A
            dy(10:14:end) = rv1.*y(4:14:end)+p.gamma*(y(8:14:end) + y(9:14:end))-(p.mu+p.rho+rv2).*y(10:14:end); % V^R_1d
            
            
            dy(11:14:end) = rv2.*y(7:14:end)-(FI.*(1-ocv.eta_2d(:,index_t))+p.mu).*y(11:14:end)+p.rho*y(14:14:end); % V^S_2d
            dy(12:14:end) = p.sigma.*FI.*(1-ocv.eta_2d(:,index_t)).*y(11:14:end)-(p.gamma+p.alpha+p.mu)*y(12:14:end); % V^I_2d
            dy(13:14:end) = (1-p.sigma).*FI.*(1-ocv.eta_2d(:,index_t)).*y(11:14:end)+rv2.*y(9:14:end)-(p.mu+p.gamma).*y(13:14:end); % V^R_2d
            dy(14:14:end) = rv1.*y(10:14:end)+p.gamma*(y(12:14:end)+y(13:14:end))-(p.mu+p.rho).*y(14:14:end); % V^R_1d

        end
    end

    
end