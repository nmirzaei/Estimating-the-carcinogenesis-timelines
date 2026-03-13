clc
clear

cancer_type='Thyroid';
Sex = 'Male';
data_lim = 'less50';
user4 = 0;
cohortmin = 1950;
max_pop = 480;

num = 500-480;
IDX1 = randi([1,max_pop],[1,num]);

csvName = strcat('Pars-cohorts-SimpleBirth\New data\',lower(cancer_type),'\',Sex,'\',data_lim,'\',num2str(cohortmin),'--',num2str(cohortmin+4),'_Pars_',cancer_type,'_',Sex,'_',data_lim,'.csv');
pp = csvread(csvName,1,0);

pp(max_pop+1:end,:)= pp(IDX1,:);

TTT = array2table(pp);
TTT.Properties.VariableNames(1:7) = {'muN','mu1','mu2','alpha1','alpha2','beta1','beta2'};
filename = strcat(num2str(cohortmin),'--',num2str(cohortmin+4),'_Pars_',cancer_type,'_',Sex,'__',data_lim,'.csv');
writetable(TTT,filename);