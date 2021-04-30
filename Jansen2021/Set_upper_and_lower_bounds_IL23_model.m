disp('Setting bounds on parameter priors.')

%% Retrieve parameter values from zero-edge model fit:
values=IL23_model_saved_chi2s.p(1,:);

%% Retrieve parameter values from 1-edge model fits:
for y=1:length(IL23_model_saved_chi2s.edges)
    if max(size(IL23_model_saved_chi2s.edges{y}))==1
        values(end+1,:)=IL23_model_saved_chi2s.p(y,:);
    end
end

%% Retrieve parameter values from full model fit:
for y=1:length(IL23_model_saved_chi2s.edges)
    if max(size(IL23_model_saved_chi2s.edges{y}))==length(edge_labels)
        values(end+1,:)=IL23_model_saved_chi2s.p(y,:);
    end
end

%% Set priors:
disp('We constrain the upper and lower bounds of our intial guess intervals.')
disp('For each parameter, we find its maximum (minimum) value from the')
disp('collections of fitted parameter values for the 0-, 1- and full-edge')
disp('models. Here, the full-edge model is the model configuration with')
disp('all signficiant edges. We set upper (lower) bounds at 2x more (10x less)')
disp('than the found maximum (minimum).')
ubs=max(values)+log10(3);
lbs=min([min(values(1:end-1,1:5)),max(values(1:end-1,6:end))],values(end,:))-1;
ar.ub(1:25)=ubs(1:25);
ar.lb(1:25)=lbs(1:25);
ar.p(1:25)=lbs(1:25)+0.1;
IL23_model_saved_chi2s.lb=ar.lb;
IL23_model_saved_chi2s.ub=ar.ub;
IL23_model_saved_chi2s.ar=ar;
arPrint
save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')