function [ocv, cati_sum] = generate_scenario(scenario_ocv, scenario_npi, time_model)

load ../data/geodata POPnodes WS_dept
POPnodes = POPnodes' * WS_dept;
POP = POPnodes';

% MAKE SCENARIO
    eff_months = 0:6:60;
    eff_days = eff_months*365.25/12;
    eff_days = round(eff_days);
    eff_adults_1d=[76, 72, 68, 0, 0, 0, 0, 0, 0, 0, 0]/100;
    eff_adults_2d=[76, 72, 68, 63, 58, 52, 46, 39, 32, 24, 15]/100;
    eff_childr_1d=[36, 34, 32, 0, 0, 0, 0, 0, 0, 0, 0]/100;
    eff_childr_2d=[36, 34, 32, 30, 27, 24, 22, 18, 15, 11, 7]/100;

    eff_1d = 0.89*eff_adults_1d+0.11*eff_childr_1d;
    eff_2d = 0.89*eff_adults_2d+0.11*eff_childr_2d;

    load '../data/ocv.mat' datespace rv_1d rv_2d eta_1d eta_2d

    %extract ocv data
    index_ocv=find(datespace==time_model(1)):find(datespace==time_model(end));
    ocv.rv_1d=rv_1d(index_ocv,:);
    ocv.rv_2d=rv_2d(index_ocv,:);
    ocv.eta_1d=eta_1d(index_ocv,:);
    ocv.eta_2d=eta_2d(index_ocv,:);


% 1 SCENARIO 0 -> NO VACCINATION
% 2 SCENARIO B -> WHAT ACTUALLY HAPPENED
% 3 SCENARIO 1 -> TWO DEPARTMENTS (CT, AT) 3.71 M DOSES OVER 730 DAYS AU PRORATA
% 4 SCENARIO 2 -> THREE DEPTS (CT, AT, OU) 9.76 M DOSES OVER 730 DAYS AU PRORATA
% 5 SCENARIO 3 -> SLOW NATIONAL 16.37 M DOSES OVER 1825 DAYS AU PRORATA
% 6 SCENARIO 4 -> FAST NATIONAL 16.37 M DOSES OVER 730 DAYS AU PRORATA
% 7 SCENARIO 5 -> HIGH COVERAGE 20.91 M DOSES OVER 730 DAYS AU PRORATA
% 8 SCENARIO b1-> TWICE AS MANY VACCINATIONS
% 9 SCENARIO b2-> VACCINATIONS STARTING ON NOV 1st

% SECOND DOSE ALWAYS GIVEN TWO WEEK AFTER THE FIRST DOSE

% SCENARIOS 1-4: 70% COVERAGE WITH TWO DOSES; 10% COVERAGE WITH ONE DOSE ->
% OF ALL DOSES, 80/150 ARE FIRST DOSES; 70/150 ARE SECOND DOSES

% SCENARIO 5: 95% COVERAGE WITH TWO DOSES; 1.67% COVERAGE WITH ONE DOSE ->
% OF ALL DOSES, 96.67/191.67 ARE FIRST DOSES; 95/191.67 ARE SECOND DOSES

% COMPUTATION OF THE EFFICIENCY AVERAGE TIME SINCE VACCINATION

switch scenario_ocv
    case 1
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        
    case 2
        ocv.rv_1d = ocv.rv_1d;
        ocv.rv_2d = ocv.rv_2d;
        ocv.eta_1d = ocv.eta_1d;
        ocv.eta_2d = ocv.eta_2d;
        
   case 8
        ocv.rv_1d = 2*ocv.rv_1d;
        ocv.rv_2d = 2*ocv.rv_2d;
        ocv.eta_1d = ocv.eta_1d;
        ocv.eta_2d = ocv.eta_2d;
        
    case 9
        load '../data/ocv.mat' datespace rv_1d rv_2d eta_1d eta_2d
        date_diff = datenum(2016,11,08) - datenum(2011,11,01) + 1;
        ocv.rv_1d = rv_1d(date_diff:date_diff+length(index_ocv)-1,:);
        ocv.rv_2d = rv_2d(date_diff:date_diff+length(index_ocv)-1,:);
        ocv.eta_1d = eta_1d(date_diff:date_diff+length(index_ocv)-1,:);
        ocv.eta_2d = eta_2d(date_diff:date_diff+length(index_ocv)-1,:);
        
        
        
    case 4
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        
        for i = 1:length(time_model)
            if (time_model(i) >= datenum(2011,10,01)) && (time_model(i) <= datenum(2011,10,01)+730)
                ocv.rv_1d(i,1) = (80/150) * 9.76e6 * POP(1)/(POP(1)+POP(2)+POP(4)) / 730;
                ocv.rv_1d(i,2) = (80/150) * 9.76e6 * POP(2)/(POP(1)+POP(2)+POP(4)) / 730;
                ocv.rv_1d(i,4) = (80/150) * 9.76e6 * POP(4)/(POP(1)+POP(2)+POP(4)) / 730;
                ocv.eta_1d(i,[1 2 4]) = max(eff_1d);
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_1d(i,[1 2 4]) = interp1(eff_days, eff_1d, time_model(i)-datenum(2011,10,01) + 730, 'linear', 'extrap');
            end
            if (time_model(i) >= datenum(2011,10,01) + 14) && (time_model(i) <= datenum(2011,10,01)+730)
                ocv.rv_2d(i,1) = (70/150) * 9.76e6 * POP(1)/(POP(1)+POP(2)+POP(4)) / (730-14);
                ocv.rv_2d(i,2) = (70/150) * 9.76e6 * POP(2)/(POP(1)+POP(2)+POP(4)) / (730-14);
                ocv.rv_2d(i,4) = (70/150) * 9.76e6 * POP(4)/(POP(1)+POP(2)+POP(4)) / (730-14);
                ocv.eta_2d(i,[1 2 4]) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) + 14)/2, 'linear', 'extrap');
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_2d(i,[1 2 4]) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) + 14)/2, 'linear', 'extrap');
            end
        end
        
    case 3
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        
        for i = 1:length(time_model)
            if (time_model(i) >= datenum(2011,10,01)) && (time_model(i) <= datenum(2011,10,01)+730)
                ocv.rv_1d(i,1) = (80/150) * 3.71e6 * POP(1)/(POP(1)+POP(2)) / 730;
                ocv.rv_1d(i,2) = (80/150) * 3.71e6 * POP(2)/(POP(1)+POP(2)) / 730;
                ocv.eta_1d(i,1:2) = max(eff_1d);
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_1d(i,1:2) = interp1(eff_days, eff_1d, time_model(i)-datenum(2011,10,01) - 730, 'linear', 'extrap');
            end
            if (time_model(i) >= datenum(2011,10,01) + 14) && (time_model(i) <= datenum(2011,10,01)+730)
                ocv.rv_2d(i,1) = (70/150) * 3.71e6 * POP(1)/(POP(1)+POP(2)) / (730-14);
                ocv.rv_2d(i,2) = (70/150) * 3.71e6 * POP(2)/(POP(1)+POP(2)) / (730-14);
                ocv.eta_2d(i,1:2) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_2d(i,1:2) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
            end
        end
                
    case 5
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        for i = 1:length(time_model)
            if (time_model(i) >= datenum(2011,10,01)) && (time_model(i) <= datenum(2011,10,01)+1825)
                for j = 1:10
                    ocv.rv_1d(i,j) = (80/150) * 16.37e6 * POP(j)/sum(POP) / 1825;
                    ocv.eta_1d(i,j) = max(eff_1d);
                end
            end
            if time_model(i) > datenum(2011,10,01) + 1825
                ocv.eta_1d(i,:) = interp1(eff_days, eff_1d, time_model(i)-datenum(2011,10,01) - 1825, 'linear', 'extrap');
            end
            if (time_model(i) >= datenum(2011,10,01) + 14) && (time_model(i) <= datenum(2011,10,01)+1825)
                for j = 1:10
                    ocv.rv_2d(i,j) = (70/150) * 16.37e6 * POP(j)/sum(POP) / (1825-14);
                    ocv.eta_2d(i,j) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
                end
            end
            if time_model(i) > datenum(2011,10,01) + 1825
                ocv.eta_2d(i,:) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
            end
        end
        
    case 6
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        for i = 1:length(time_model)
            if (time_model(i) >= datenum(2011,10,01)) && (time_model(i) <= datenum(2011,10,01)+730)
                for j = 1:10
                    ocv.rv_1d(i,j) = (80/150) * 16.37e6 * POP(j)/sum(POP) / 730;
                    ocv.eta_1d(i,j) = max(eff_1d);
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_1d(i,:) = interp1(eff_days, eff_1d, time_model(i)-datenum(2011,10,01) - 730, 'linear', 'extrap');
            end
            if (time_model(i) >= datenum(2011,10,01) + 14) && (time_model(i) <= datenum(2011,10,01)+730)
                for j = 1:10
                    ocv.rv_2d(i,j) = (70/150) * 16.37e6 * POP(j)/sum(POP) / (730-14);
                    ocv.eta_2d(i,j) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_2d(i,:) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
            end
        end
        
    case 7
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        for i = 1:length(time_model)
            if (time_model(i) >= datenum(2011,10,01)) && (time_model(i) <= datenum(2011,10,01)+730)
                for j = 1:10
                    ocv.rv_1d(i,j) = (96.67/191.67) * 20.91e6 * POP(j)/sum(POP) / 730;
                    ocv.eta_1d(i,j) = max(eff_1d);
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_1d(i,:) = interp1(eff_days, eff_1d, time_model(i)-datenum(2011,10,01) - 730, 'linear', 'extrap');
            end
            if (time_model(i) >= datenum(2011,10,01) + 14) && (time_model(i) <= datenum(2011,10,01)+730)
                for j = 1:10
                    ocv.rv_2d(i,j) = (95/191.67) * 20.91e6 * POP(j)/sum(POP) / (730-14);
                    ocv.eta_2d(i,j) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_2d(i,:) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
            end
        end
        
        case 10
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        for i = 1:length(time_model)
            if (time_model(i) >= datenum(2011,10,01)) && (time_model(i) <= datenum(2011,10,01)+730)
                for j = 1:10
                    ocv.rv_1d(i,j) = 0.95* POP(j) / 730;
                    ocv.eta_1d(i,j) = max(eff_1d);
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_1d(i,:) = interp1(eff_days, eff_1d, time_model(i)-datenum(2011,10,01) - 730, 'linear', 'extrap');
            end
            if (time_model(i) >= datenum(2011,10,01) + 14) && (time_model(i) <= datenum(2011,10,01)+730 + 14)
                for j = 1:10
                    ocv.rv_2d(i,j) = 0.95 * POP(j) / (730);
                    ocv.eta_2d(i,j) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730 + 14
                ocv.eta_2d(i,:) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
            end
        end
        
        case 11
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        for i = 1:length(time_model)
            if (time_model(i) >= datenum(2011,10,01)) && (time_model(i) <= datenum(2011,10,01)+730)
                for j = 1:10
                    ocv.rv_1d(i,j) = 0.90* POP(j) / 730;
                    ocv.eta_1d(i,j) = max(eff_1d);
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_1d(i,:) = interp1(eff_days, eff_1d, time_model(i)-datenum(2011,10,01) - 730, 'linear', 'extrap');
            end
            if (time_model(i) >= datenum(2011,10,01) + 14) && (time_model(i) <= datenum(2011,10,01)+730 + 14)
                for j = 1:10
                    ocv.rv_2d(i,j) = 0.90 * POP(j) / (730);
                    ocv.eta_2d(i,j) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730 + 14
                ocv.eta_2d(i,:) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
            end
        end
        
        case 12
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        for i = 1:length(time_model)
            if (time_model(i) >= datenum(2011,10,01)) && (time_model(i) <= datenum(2011,10,01)+730)
                for j = 1:10
                    ocv.rv_1d(i,j) = 0.75* POP(j) / 730;
                    ocv.eta_1d(i,j) = max(eff_1d);
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_1d(i,:) = interp1(eff_days, eff_1d, time_model(i)-datenum(2011,10,01) - 730, 'linear', 'extrap');
            end
            if (time_model(i) >= datenum(2011,10,01) + 14) && (time_model(i) <= datenum(2011,10,01)+730 + 14)
                for j = 1:10
                    ocv.rv_2d(i,j) = 0.75 * POP(j) / (730);
                    ocv.eta_2d(i,j) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730 + 14
                ocv.eta_2d(i,:) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
            end
        end
        
        case 13
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        for i = 1:length(time_model)
            if (time_model(i) >= datenum(2011,10,01)) && (time_model(i) <= datenum(2011,10,01)+730)
                for j = 1:10
                    ocv.rv_1d(i,j) = 0.50* POP(j) / 730;
                    ocv.eta_1d(i,j) = max(eff_1d);
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_1d(i,:) = interp1(eff_days, eff_1d, time_model(i)-datenum(2011,10,01) - 730, 'linear', 'extrap');
            end
            if (time_model(i) >= datenum(2011,10,01) + 14) && (time_model(i) <= datenum(2011,10,01)+730 + 14)
                for j = 1:10
                    ocv.rv_2d(i,j) = 0.50 * POP(j) / (730);
                    ocv.eta_2d(i,j) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730 + 14
                ocv.eta_2d(i,:) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
            end
        end
        
        case 14
        ocv.rv_1d = zeros(length(time_model),10);
        ocv.rv_2d = zeros(length(time_model),10);
        ocv.eta_1d = zeros(length(time_model),10);
        ocv.eta_2d = zeros(length(time_model),10);
        for i = 1:length(time_model)
            if (time_model(i) >= datenum(2011,10,01)) && (time_model(i) <= datenum(2011,10,01)+730)
                for j = 1:10
                    ocv.rv_1d(i,j) = 0.20* POP(j) / 730;
                    ocv.eta_1d(i,j) = max(eff_1d);
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730
                ocv.eta_1d(i,:) = interp1(eff_days, eff_1d, time_model(i)-datenum(2011,10,01) - 730, 'linear', 'extrap');
            end
            if (time_model(i) >= datenum(2011,10,01) + 14) && (time_model(i) <= datenum(2011,10,01)+730 + 14)
                for j = 1:10
                    ocv.rv_2d(i,j) = 0.20 * POP(j) / (730);
                    ocv.eta_2d(i,j) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
                end
            end
            if time_model(i) > datenum(2011,10,01) + 730 + 14
                ocv.eta_2d(i,:) = interp1(eff_days, eff_2d, (time_model(i)-datenum(2011,10,01) - 14)/2, 'linear', 'extrap');
            end
        end
end

    % transpose
    ocv.rv_1d=ocv.rv_1d';
    ocv.rv_2d=ocv.rv_2d';
    ocv.eta_1d=ocv.eta_1d';
    ocv.eta_2d=ocv.eta_2d';

% SCENARIO 0 -> NO NPI
% SCENARIO B -> WHAT ACTUALLY HAPPENED
% SCENARIO 1 -> TWICE AS NPI
% SCENARIO 2 -> TENFOLD AS NPI
% SCENARIO 3 -> 1 NPI EVERY 100 CASES
% SCENARIO 4 -> 1 NPI EVERY 50 CASES
% SCENARIO 5 -> 1 NPI EVERY 10 CASES
% SCENARIO 1b -> ANTICIPATED NPI
% SCENARIO 1.5 -> 5FOLD AS NPI

%extract cati data
    if scenario_npi == 2 || scenario_npi == 3 ||  scenario_npi == 4 || scenario_npi == 8 || scenario_npi == 9
        load ../data/wash date_list_cati cati_day
        date_list_cati = double(date_list_cati)';
        index_cati=find(date_list_cati==time_model(1)):find(date_list_cati==time_model(end));
        cati_sum = cumsum(cati_day,1);
        cati_sum = double(cati_sum(index_cati,:));
    end

    cases_week = csvread('../data/cases.csv',1,1)'; 

switch scenario_npi
    case 1
        cati_sum = zeros(length(time_model),10);
    case 2
        cati_sum = cati_sum;
    case 3
        cati_sum = 2*cati_sum;
    case 4
        cati_sum = 10*cati_sum;
    case 9
        cati_sum = 5*cati_sum;
    case 5
        cati_day = zeros(length(time_model),10);
        for k = 1:size(cases_week,2)
            time_base = datenum(2010,10,23) + (k-1)*7;
            datestr(time_base)
            if time_base >= datenum(2011,10,01)
                for dep = 1:10
                    cati_day(3+7*(k-1),dep) = round(cases_week(dep,k-1)*0.01);
                end
            end
        end
        cati_day(isnan(cati_day)) = 0;
        cati_sum = cumsum(cati_day,1);
        
    case 6
        cati_day = zeros(length(time_model),10);
        for k = 1:size(cases_week,2)
            time_base = datenum(2010,10,23) + (k-1)*7;
            datestr(time_base)
            if time_base >= datenum(2011,10,01)
                for dep = 1:10
                    cati_day(3+7*(k-1),dep) = round(cases_week(dep,k-1)*0.02);
                end
            end
        end
        cati_day(isnan(cati_day)) = 0;
        cati_sum = cumsum(cati_day,1);
        
    case 7
        cati_day = zeros(length(time_model),10);
        for k = 1:size(cases_week,2)
            time_base = datenum(2010,10,23) + (k-1)*7;
            datestr(time_base)
            if time_base >= datenum(2011,10,01)
                for dep = 1:10
                    cati_day(3+7*(k-1),dep) = round(cases_week(dep,k-1)*0.1);
                end
            end
        end
        cati_day(isnan(cati_day)) = 0;
        cati_sum = cumsum(cati_day,1);
        
    case 8
        cati_sum2 = cati_sum;
        cati_sum = zeros(length(time_model),10);
        cati_sum(346:1802,:) = cati_sum2(991:end,:);
        for i = 1803:size(cati_sum,1);
            cati_sum(i,:) = cati_sum(1802,:);
        end
        
end
end
