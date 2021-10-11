function [ModPred,times_out, y] = m4c(x,Extra)

    load data/precipitation rainfall_day date_list
    load data/geodata nnodes POPnodes WS_dept dist_road

    %Fitting Parameters
    p.theta= x(1);
    p.m= x(2);
    p.D= x(3);
    p.phi= x(4);
    p.rho = x(5);
    p.sigma= x(6);
    p.muB= x(7);
    p.beta0= x(8);
    p.psi = x(9);
    p.t0= round(x(10));    
    p.r = x(11);

    %pre-defined parameters
    p.gamma=0.2;               %rate at which people recover from cholera (day^-1)
    p.mu=1/(61.4*365);         %population natality and mortality rate (day^-1)
    p.alpha=-log(0.98)*p.gamma;  %mortality rate due to cholera (day^-1)
    
    %set time
    t0=0;                    %delay between the first day of data available and the day of the onset of the epidemic (initial condition)
    t_initial=datenum('20.10.2010','dd.mm.yyyy');
    t_final=datenum('30.06.2013','dd.mm.yyyy');
    time_data=t_initial+(7-weekday(t_initial))+(0:141-1)*7; 
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
    
    % initial condition
    nnodes = 10;
    POPnodes = (POPnodes' * WS_dept)';
    
    initial_points=[1 2];
    I_initial=[1100 1000]/p.sigma;
    
    y0=zeros(6,nnodes);
    y0(1,:)=POPnodes; %whole population
    y0(1,initial_points) = y0(1,initial_points)-I_initial; %susceptible = pop minus first infected
    y0(2,initial_points) = p.sigma*I_initial; %first symptomatic infections
    y0(3,initial_points) = (1-p.sigma)*I_initial; %first removed
    y0(6,initial_points) = p.sigma*I_initial; %reported cases @ beginning
    y0(5,initial_points) = (y0(2,initial_points)+p.r*y0(3,initial_points))*p.theta/p.muB./POPnodes(initial_points)';

    %run model
    y=SIB(p,nnodes,POPnodes,fluxes,rainfall_day,time_model,y0);
    cumcases_AD1_t=y(:,6:6:end)';

    %upscale from daily to weekly time scale
    index_time=zeros(length(time_data),1);
    for iii=1:length(time_data)
        index_time(iii)=find(time_model==time_data(iii));
    end

    cumcases_AD1_week=cumcases_AD1_t(:,index_time);
    cases_AD1_week=diff([zeros(size(cumcases_AD1_week,1),1) cumcases_AD1_week],1,2); %first column is 0, makes the difference between each element along the row.
    
    ModPred = cases_AD1_week(:,1:141);
    times_out = time_data(1:141);
    %%% ODE PART
    
    function ytop=SIB(p,nnodes,H,fluxes,rainfall_day,tspan,y0)

        step_size = 7;
        for iii = 1:step_size:length(tspan)
            try
                tspan_sub = tspan(iii:iii+step_size);
            catch
                clear tspan_sub
                tspan_sub = tspan(iii:end);
            end
        
        
            [~,y]=ode45(@eqs,tspan_sub,y0);

            try
                ytop(iii:iii+step_size,:) = y;
            catch
                ytop(iii:length(tspan),:) = y;
            end
            y0 = reshape(y(end,:),6,nnodes);
            
            clear y
        end

        function dy=eqs(t,y)
            index_t=floor(t-tspan(1))+1;

            dy=zeros(6*nnodes,1);

            temp=fluxes*(y(5:6:end)./(y(5:6:end)+1));    
            

            if index_t<p.t0+1
                Ceff = y(6:6:end);
            else
                Ceff = y(6:6:end)-ytop(index_t-p.t0,6:6:end)';
            end
            beta_t = p.beta0*exp(-Ceff./H/p.psi);
            
            theta_t = p.theta*(1+p.phi*rainfall_day(:,index_t));

            %TT = FI*S
            TT=((1-p.m)*y(5:6:end)./(1+y(5:6:end))+temp*p.m).*(beta_t).*y(1:6:end);
            GI = y(2:6:end)+p.r*((1-p.m)*y(3:6:end) + p.m*fluxes'*y(3:6:end));
            
            dy(1:6:end)=p.mu*(H-y(1:6:end))- TT +p.rho*y(4:6:end);        %S
            dy(2:6:end)=p.sigma.*TT - (p.gamma+p.alpha+p.mu)*y(2:6:end);       %I
            dy(3:6:end)=(1-p.sigma).*TT - (p.gamma+p.mu)*y(3:6:end);       %I
            dy(4:6:end)=p.gamma*(y(2:6:end)+y(3:6:end)) - (p.mu+p.rho).*y(4:6:end);       %R
            dy(5:6:end)=-(p.muB).*y(5:6:end)+theta_t.*GI./H;  %B
            dy(6:6:end)=p.sigma.*TT;   %Cum(I)  

        end
    end

    
end