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
cancer_type='CRC';
Sex = input('Please type in the gender (Male, Female):','s');
data_lim = input('Enter age cutoff (less50, less70):','s');
year = input('Enter the desired 5-year cohort lower bound:');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Doubling rates and stem cells number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bm = 2.572034379774019;
N0 = 2E8; 
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
% Solve for best-fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TT = min(T,max(Age));
tt_plot = 0:min(T,TT);
x_ic_best = [1 0 1 0 1 0 1 0];
[~,yp_best] = ode15s(@(t,x) hazardfunc_multi_malignant_cells_SimpleBirth(t,x,pp(1,:),Cell_num_array(:,3),N0,bm), ...
                     tt_plot, x_ic_best);
y_best = yp_best(:,2);  % malignant hazard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
IDX = find(Age==max(tt_plot));
plot(tt_plot, y_best, 'k--', 'LineWidth', 1.5,'Color',[0.6 0.8 0.8]);
plot(Age(1:IDX),data1(1:IDX),'.','MarkerSize',12,'Color','black')
xlabel('Age (years)')
ylabel('Malignant hazard rate')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit Diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outIC= aicbic_msce_poisson(pp(1,:), T, Age, data1, Cell_num_array(:,3), N0, bm);
fprintf('AIC:')
disp(outIC.AIC)
fprintf('BIC:')
disp(outIC.BIC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%