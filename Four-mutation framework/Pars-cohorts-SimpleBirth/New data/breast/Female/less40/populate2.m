clc
clear

cancer_type='Breast';
Sex = 'Female';
data_lim = 'less40';
user4 = 0;
cohortmin = 1965;
max_pop = 500;

num = 600-500;
IDX1 = randi([1,500],[1,num]);

csvName = strcat('Pars-cohorts-SimpleBirth\New data\',lower(cancer_type),'\',Sex,'\',data_lim,'\',num2str(cohortmin),'--',num2str(cohortmin+4),'_Pars_',cancer_type,'_',Sex,'_',data_lim,'.csv');
pp = csvread(csvName,1,0);

pp(IDX1,:)= pp(501:600,:);
pp(500+1:end,:)=[];
TTT = array2table(pp);
TTT.Properties.VariableNames(1:10) = {'muN','mu1','mu2','mu3','alpha1','alpha2','alpha3','beta1','beta2','beta3'};
filename = strcat(num2str(cohortmin),'--',num2str(cohortmin+4),'_Pars_',cancer_type,'_',Sex,'_',data_lim,'.csv');
writetable(TTT,filename);