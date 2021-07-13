%% Load model:
arLoadModel('IL23_model');

%% Load data:
arLoadData('Data_IL23_model')

%% Compile:
arCompileAll();

%% Set configuration options:
ar.config.fiterrors = -1; % no error fitting
ar.config.ploterrors = -2;
ar.config.useLHS = true;
ar.config.restartLHS = true;

%% Set parameters:
arSetPars(ar.pLabel,0,1,1,-20,10);
arSetPars('sd_TNFA_mfi',0,2,[],0,0);
arSetPars('sd_IL1A_mfi',0,2,[],0,0);
arSetPars('sd_IL1B_mfi',0,2,[],0,0);
arSetPars('sd_IL6_mfi',0,2,[],0,0);
arSetPars('sd_IL10_mfi',0,2,[],0,0);
arSetPars('sd_p19_mfi',0,2,[],0,0);

%% Compute variance from data:
data_var = [];
for i=1:size(vertcat(ar.model.data.name),1)
    data_var=[data_var;ar.model.data(i).yExp-mean(ar.model.data(i).yExp,'omitnan')];
end
data_var=sqrt(var(data_var,'omitnan'));

for j=1:size(vertcat(ar.model.data.name),1)
    ar.model.data(j).yExpStd=repmat(data_var,size(ar.model.data(j).yExpStd,1),1);
end

%% Save:
save('IL23_model.mat','ar');