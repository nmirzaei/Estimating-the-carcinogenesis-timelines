clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Turning off warnings.
%Adding the functions and data path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off', 'all')
addpath("Data\")
addpath("Functions\")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%User input for cancer type and gender
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cancer_type='Breast';
Sex = 'Female';
data_lim = input('Enter age cutoff (less40, less50, less70):','s');
year = input('Enter the desired 5-year cohort lower bound:');
Menarche_Age = 12.69;
k_lim=500; %Number of estimated parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Doubling rate and number of stem cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bm = 5.244100617548181;
N0 = 1.74E10; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(data_lim,'less40')
    T=39;
elseif strcmp(data_lim,'less50')
    T=49;
elseif strcmp(data_lim,'less70')
    T=70;
else
    error('The age cutoff has to be entered as less40, less50 or less70')
end
t = 0:T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reading data and sorting by age
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat('Data\New data\',lower(cancer_type),'\data1_5Yr_Age_',Sex,'_mixed_seer.csv');
data = csvread(filename,1,0);

AgeSortData = sortrows(data,3);  %sort all rows based on the age column 
cohortmin = min(AgeSortData(:,1)); %Min cohort year
cohortmax = max(AgeSortData(:,1)); %Max cohort year
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tumor Size data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = strcat('Pars-cohorts-SimpleBirth\New data\',lower(cancer_type),'\size.csv');
Size_data = csvread(name,1,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following lines are to extract the age and cell sizes for the 
%current cohort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cohortmin=year; 
idx = find(AgeSortData(:,1)==cohortmin);
data1 = AgeSortData(idx,6);
Yeardx = AgeSortData(idx,4);
Age = AgeSortData(idx,3);
Agemin = min(Age);
idx1 = find(Size_data(:,1)==Yeardx(1,:));
idx2 = find(Size_data(:,1)==Yeardx(end,:));
Cell_num_array =[Size_data(idx1-Agemin:idx2,1) Size_data(idx1-Agemin:idx2,2) Size_data(idx1-Agemin:idx2,3)];    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameter read
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
csvName = strcat('Pars-cohorts-SimpleBirth\New data\',lower(cancer_type),'\',Sex,'\',data_lim,'\',num2str(cohortmin),'--',num2str(cohortmin+4),'_Pars_',cancer_type,'_',Sex,'_',data_lim,'.csv');
pp = csvread(csvName,1,0);
PARAMS = pp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the timelines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:k_lim
    TT = min(T,max(Age));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Model (A)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_ic1 = [1 PARAMS(k,1)*N0];
    [tp1,yp1] = ode15s(@(t,x,p) hazardfunc_sub_model_1_brc(t,x,p,N0,Menarche_Age),0:TT, x_ic1, [], PARAMS(k,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Model (B)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_ic2 = [1 0 1 -PARAMS(k,2)];
    [tp2,yp2] = ode15s(@(t,x,p) hazardfunc_sub_model_2_brc(t,x,p,N0,Menarche_Age),0:TT, x_ic2, [], PARAMS(k,:)); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Model (C)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_ic3 = [1 0 1 0 1 -PARAMS(k,3)];
    [tp3,yp3] = ode15s(@(t,x,p) hazardfunc_sub_model_3_brc(t,x,p,N0,Menarche_Age),0:TT, x_ic3, [], PARAMS(k,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Model (D)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_ic4 = [1 0 1 0 1 0 1 -PARAMS(k,4)];
    [tp4,yp4] = ode15s(@(t,x,p) hazardfunc_single_malignant_cells_brc(t,x,p,N0,Menarche_Age),0:TT, x_ic4, [], PARAMS(k,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extracting times after menarche
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_span  =0:min(T,TT);
    idx_delay = find(t_span >= Menarche_Age-1);
    t_delayed = t_span(idx_delay);
    dt = t_span(2) - t_span(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Timeline calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        P0_prime_1 = [0;yp1(idx_delay(2:end),2).*yp1(idx_delay(2:end),1)];
    
        P0_prime_2 = yp2(idx_delay,2).*yp2(idx_delay,1);
        P1_prime_2 = -yp2(idx_delay,4);
    
        P0_prime_3 = yp3(idx_delay,2).*yp3(idx_delay,1);
        P1_prime_3 = -yp3(idx_delay,4);
        P2_prime_3 = -yp3(idx_delay,6);
    
        P0_prime_4 = yp4(idx_delay,2).*yp4(idx_delay,1);
        P1_prime_4 = -yp4(idx_delay,4);
        P2_prime_4 = -yp4(idx_delay,6);
        P3_prime_4 = -yp4(idx_delay,8);
    
        % Normalize to get densities
        f_FSMC = P0_prime_1 / trapz(t_delayed, P0_prime_1); 
    
    
        f_SSMC_given_FSMC = P1_prime_2 / trapz(t_delayed, P1_prime_2); 
        f_SSMC = conv(f_FSMC,f_SSMC_given_FSMC)*dt; 
        f_SSMC = f_SSMC(1:length(t_delayed));
        f_SSMC = f_SSMC/trapz(t_delayed,f_SSMC);
    
    
        f_TSMC_given_SSMC = P2_prime_3 / trapz(t_delayed, P2_prime_3); 
        f_TSMC = conv(f_SSMC,f_TSMC_given_SSMC)*dt;
        f_TSMC = f_TSMC(1:length(t_delayed));
        f_TSMC = f_TSMC/trapz(t_delayed,f_TSMC);
    
        f_malignant_given_TSMC = P3_prime_4 / trapz(t_delayed, P3_prime_4);
        f_malignant = conv(f_TSMC,f_malignant_given_TSMC)*dt;
        f_malignant = f_malignant(1:length(t_delayed));
        f_malignant = f_malignant/trapz(t_delayed,f_malignant);
        
    
        % Expected ages (mean of densities)
        E_FSMC(k,1) = trapz(t_delayed, t_delayed .* f_FSMC');
        E_SSMC(k,1) = trapz(t_delayed, t_delayed.* f_SSMC');
        E_TSMC(k,1) = trapz(t_delayed, t_delayed.* f_TSMC');
        E_Malignant(k,1) = trapz(t_delayed, t_delayed .* f_malignant');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch
        % Expected ages (mean of densities)
        E_FSMC(k,1) = NaN;
        E_SSMC(k,1) = NaN;
        E_TSMC(k,1) = NaN;
        E_Malignant(k,1) = NaN;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TTT = array2table([E_FSMC E_SSMC E_TSMC E_Malignant]);
TTT.Properties.VariableNames(1:4) = {'Stem to FSMC','FSMC to SSMC','SSMC to TSMC','TSMC to Malignant'};
data = TTT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize result container
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varNames = data.Properties.VariableNames;
nCols = width(data);

fprintf('Column\t\tMean (CI: lower-upper)\n');
fprintf('----------------------------------------\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the Mean and CIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nCols
    columnData = data{:, i};
    columnData = columnData(~isnan(columnData));

    meanVal = mean(columnData);
    q = prctile(columnData, [2.5 50 97.5]); 

    fprintf('%s:\tmean %.1f | median %.1f | 95%% UI: %.1f–%.1f\n', ...
        varNames{i}, meanVal, q(2), q(1), q(3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

