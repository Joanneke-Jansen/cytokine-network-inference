%% Initialise:
disp('Initialise:')
clear all; close all;
arInit

%% Compile the model:
disp('Compile:')
try load('IL23_model.mat','ar')
catch
    Compile_IL23_model
    load('IL23_model.mat','ar')
end

%% Fit the model with all significant edges to the data (the N-sized model):
disp('Fit the model with all significant edges to the data:')
TSTART = tic;
%Elapsed time is 950.178964 seconds.
First_fit_IL23_model

%% Compute D and AIC for the initial list of configurations (Method Step 2):
Initial_list_IL23_model
[~,n_o]=min(IL23_model_saved_chi2s.AIC);
disp('The model with edges:')
disp((IL23_model_saved_chi2s.initial_model{n_o}))
disp('currently has the lowest relative AIC')

%% We check whether the initial model of size n_o is also the optimal model of size n_o (Method Step 3):
[IL23_model_saved_chi2s]=Procedure_Omega_IL23_model(n_o,IL23_model_saved_chi2s);

%% We check whether a network configuration exists, with a lower relative AIC than S_selected:
% for every model size n=1, 2, ..., n_o+p_{max}. (Method Step 4):
p_max=floor((IL23_model_saved_chi2s.chi2s(n_o)-IL23_model_saved_chi2s.chi2s(end))/2);
for I=1:(n_o+p_max)
    [IL23_model_saved_chi2s]=Procedure_Omega_hat_IL23_model(I,IL23_model_saved_chi2s);
end
toc(TSTART)

%% We conclude:
disp('The model with edges:')
[~,n_s]=min(IL23_model_saved_chi2s.AIC);
selected_configuration=IL23_model_saved_chi2s.initial_model{n_s};
disp(selected_configuration)
disp(' has the lowest relative AIC')

for y=1:length(IL23_model_saved_chi2s.edges)
    if isequal(sort(IL23_model_saved_chi2s.edges{y}),sort(IL23_model_saved_chi2s.initial_model{n_s}))
        index=y;
        ar.p=IL23_model_saved_chi2s.p(index,:);
    end
end

disp('With edge values')
for i=1:length(selected_configuration)
    disp([selected_configuration{i},' = ', num2str(10^arGetPars(selected_configuration{i}))])
end
