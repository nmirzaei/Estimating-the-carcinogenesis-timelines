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
cancer_type='Thyroid';
Sex = input('Please type in the gender (Male, Female):','s');
data_lim = input('Enter age cutoff (less50, less70):','s');
year = input('Enter the desired 5-year cohort lower bound:');
k_lim=500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Doubling rates and stem cells number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bm = 1.458497192428218;
N0 = 6.5E7; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(data_lim,'less50')
    T=49;
elseif strcmp(data_lim,'less70')
    T=70;
else
    error('The age cutoff has to be entered as less50 or less70')
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
    [tp1,yp1] = ode15s(@(t,x,p) hazardfunc_sub_model_1(t,x,p,N0),0:TT, x_ic1, [], PARAMS(k,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Model (B)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_ic2 = [1 0 1 -PARAMS(k,2)];
    [tp2,yp2] = ode15s(@(t,x,p) hazardfunc_sub_model_2(t,x,p,N0),0:TT, x_ic2, [], PARAMS(k,:));    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Model (C)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_ic3 = [1 0 1 0 1 -PARAMS(k,3)];
    [tp3,yp3] = ode15s(@(t,x,p) hazardfunc_single_malignant_cells(t,x,p,N0),0:TT, x_ic3, [], PARAMS(k,:));   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Timeline calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_span  = 0:min(T,TT);
    dt = t_span(2) - t_span(1);
    P0_prime_1 = [0;yp1(2:min(T,TT)+1,2).*yp1(2:min(T,TT)+1,1)];

    P0_prime_2 = yp2(1:min(T,TT)+1,2).*yp2(1:min(T,TT)+1,1);
    P1_prime_2 = -yp2(1:min(T,TT)+1,4);

    P0_prime_3 = yp3(1:min(T,TT)+1,2).*yp3(1:min(T,TT)+1,1);
    P1_prime_3 = -yp3(1:min(T,TT)+1,4);
    P2_prime_3 = -yp3(1:min(T,TT)+1,6);

    % Normalize to get densities
    f_FSMC = P0_prime_1 / trapz(t_span, P0_prime_1); 

    f_SSMC_given_FSMC = P1_prime_2 / trapz(t_span, P1_prime_2);
    f_SSMC = conv(f_FSMC,f_SSMC_given_FSMC)*dt;
    f_SSMC = f_SSMC(1:length(t_span));
    f_SSMC = f_SSMC/trapz(t_span,f_SSMC);


    f_Malignant_given_SSMC = P2_prime_3 / trapz(t_span, P2_prime_3); 
    f_Malignant = conv(f_SSMC,f_Malignant_given_SSMC)*dt;
    f_Malignant = f_Malignant(1:length(t_span));
    f_Malignant = f_Malignant/trapz(t_span,f_Malignant);

    % Expected ages (mean of densities)
    E_FSMC(k,1) = trapz(t_span, t_span .* f_FSMC');
    E_SSMC(k,1) = trapz(t_span, t_span .* f_SSMC');
    E_Malignant(k,1) = trapz(t_span, t_span .* f_Malignant');
    E_FSMC(E_FSMC<0) = 0;
    E_SSMC(E_SSMC<0) = 0;
    E_Malignant(E_Malignant<0) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TTT = array2table([E_FSMC E_SSMC E_Malignant]);
TTT.Properties.VariableNames(1:3) = {'Stem to FSMC','FSMC to SSMC','SSMC to Malignant'};
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
for i = 1:nCols
    columnData = data{:, i};
    columnData = columnData(~isnan(columnData));

    meanVal = mean(columnData);
    q = prctile(columnData, [2.5 50 97.5]);   

    fprintf('%s:\tmean %.1f | median %.1f | 95%% UI: %.1f–%.1f\n', ...
        varNames{i}, meanVal, q(2), q(1), q(3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%