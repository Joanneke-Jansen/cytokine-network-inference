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

%% Fit the 0-edge, 1-edge, and model with all significant edges to the data:
disp('Fit the 0-edge, 1-edge, and model with all significant edges to the data:')
tic;
First_fits_IL23_model
first_time=toc;

%% Set upper and lower bounds on the intervals from which parameter initial guesses are sampled:
Set_upper_and_lower_bounds_IL23_model

%% Create initial list of configurations:
tic;
Initial_list_IL23_model
second_time=toc;

%% Compute AIC for initial list of configurations:
for w=1:length(edge_labels)
    index=0;
    min_chi2=IL23_model_saved_chi2s.chi2(1);
    for y=1:length(IL23_model_saved_chi2s.edges)
        if size(IL23_model_saved_chi2s.edges{y},2)==w
            chi2=IL23_model_saved_chi2s.chi2(y);
            if chi2<min_chi2
                index=y;
                min_chi2=chi2;
            end
        end
    end
    IL23_model_saved_chi2s.initial_model{w}=IL23_model_saved_chi2s.edges{index};
    IL23_model_saved_chi2s.chi2s(w)=min_chi2;
end

for i=1:size(IL23_model_saved_chi2s.initial_model,2)
    number_of_edges=size(IL23_model_saved_chi2s.initial_model{i},2);
    IL23_model_saved_chi2s.AIC(i)=IL23_model_saved_chi2s.chi2s(i)+2*number_of_edges;
end

IL23_model_saved_chi2s.AIC=IL23_model_saved_chi2s.AIC-min(IL23_model_saved_chi2s.AIC);
[~,n_o]=min(IL23_model_saved_chi2s.AIC);

disp('The model with edges:')
IL23_model_saved_chi2s.initial_model{n_o}
disp('currently has the lowest relative AIC')

%% We check whether the initial model of size n_o is also the optimal model of size n_o:
[IL23_model_saved_chi2s]=Check_all_configurations_IL23_model(n_o,IL23_model_saved_chi2s);

%% If no better model is found, we continue:
p_max=floor((IL23_model_saved_chi2s.chi2s(n_o)-IL23_model_saved_chi2s.chi2s(end))/2);
disp(['We can add ', num2str(p_max), ' edges:'])
for i=1:p_max
    [IL23_model_saved_chi2s]=Check_all_configurations_IL23_model(n_o+p_max,IL23_model_saved_chi2s);
end

%% If no better model is found, we continue:
n_omega=n_o;
disp(['We substract 1 to psi_max edges:'])
psi=0;
while IL23_model_saved_chi2s.chi2s(n_omega-psi)-IL23_model_saved_chi2s.chi2s(n_omega)<2*n_omega
     psi=psi+1;
    [IL23_model_saved_chi2s]=Check_all_configurations_IL23_model(n_omega-psi,IL23_model_saved_chi2s);
    if psi==n_omega
        return
    end
end
psi_max=psi;
disp(['We substracted up to ', num2str(psi_max), ' edges.'])

%% If no better model is found, we conclude:
disp('The model with edges:')
selected_configuration=IL23_model_saved_chi2s.initial_model{n_omega};
disp(selected_configuration)
disp(' has the lowest relative AIC')

for y=1:length(IL23_model_saved_chi2s.edges)
        if isequal(sort(IL23_model_saved_chi2s.edges{y}),sort(IL23_model_saved_chi2s.initial_model{n_omega}))
                index=y;
                ar.p=IL23_model_saved_chi2s.p(index,:);
        end
end

disp('With edge values')
for i=1:length(selected_configuration)
    disp([selected_configuration{i},' = ', num2str(10^arGetPars(selected_configuration{i}))])
end
