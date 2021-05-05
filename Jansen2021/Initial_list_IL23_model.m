%% We calculate D for the initial list of network configurations

disp('Retrieve parameter values of N-sized model from struct:')
for y=1:length(IL23_model_saved_chi2s.edges)
    if max(size(IL23_model_saved_chi2s.edges{y}))==max(size(edge_labels))
        ar.p=IL23_model_saved_chi2s.p(y,:)';
        flag=1;
    end
end

disp('Sort parameter values in descending order:')
edge_lables_sorted=ar.pLabel(contains(ar.pLabel,edge_labels));
[~,edgevalues_sorted]=sort(ar.p(contains(ar.pLabel,edge_labels)),'descend');
edge_lables_sorted=edge_lables_sorted(edgevalues_sorted);

disp('Compute D for the initial list of network configurations:')
for n=1:length(edge_labels)
    
    %The initial n-sized model consists of the n edges with the
    %largest fitted values of the N-sized model.
    subconfiguration=edge_lables_sorted(1:n);
    disp('Model with edges:')
    disp(subconfiguration)
    
    %Test if already computed:
    flag=0;
    for y=1:length(IL23_model_saved_chi2s.edges)
        if size(IL23_model_saved_chi2s.edges{y},2)==size(subconfiguration,2)
            if sum(contains(IL23_model_saved_chi2s.edges{y},subconfiguration))==size(subconfiguration,2)
                disp('Chi2 already computed.')
                flag=1;
            end
        end
    end
    
    if ~flag
    arSetPars('LPS_IL10_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL1A_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL1B_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL23_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL6_c',0,1,[],-1.5,0);
    arSetPars('LPS_TNFA_c',0,1,[],-1.5,0);
    
        %Remove all edges:
        for p=1:length(edge_labels)
            arSetPars(edge_labels{p},-100,2,[],-100,[]);
        end
        %Insert edges of interest:
        for z=1:size(subconfiguration,2)
            arSetPars(subconfiguration{z},-1,1,[],-2,1);
        end
        
        arFitLHS_Jansen2021(51);
        
        if sum([((ar.p(ar.qFit==1)-ar.lb(ar.qFit==1))<=0.01),((ar.ub(ar.qFit==1)-ar.p(ar.qFit==1))<=0.01)])
            %disp('Warning: returned a parameter value outside of initial-guess sampling interval.')
        end
        
        if ar.LhsSampleSizeCalculation.D > 1
            IL23_model_saved_chi2s.multiple_minima_number_of_minima(end+1)=ar.LhsSampleSizeCalculation.D;
            IL23_model_saved_chi2s.multiple_minima_number_of_minima_observed_once(end+1)=ar.LhsSampleSizeCalculation.f1;
            IL23_model_saved_chi2s.multiple_minima_edges{end+1}=subconfiguration;
            if max(min(ar.ps(ar.chi2s-min(ar.chi2s)>1e-4,ar.qFit==1)'))>-16
                %Comment out the following error to increase the number
                %of initial guesses until all local minima are observed more
                %than once (such that the estimated number of minima by the chao2 method matches the observed number of minima):
                error('More than one minimum found with all parameter values larger than 10^(-16).')
                if ar.LhsSampleSizeCalculation.f1>0
                    temp=arLhsSampleSizeCalculation_Jansen2021;
                    chi2s=ar.chi2s; %Save old values
                    while temp.f1>0
                        arFitLHS_Jansen2021(51); %Compute 51 more
                        ar.chi2s=[ar.chi2s,chi2s]; %Combine
                        temp=arLhsSampleSizeCalculation_Jansen2021;
                        chi2s=ar.chi2s; %Save old values
                        ar.chi2s=ar.chi2s(1:51);
                    end
                end
            end
        end
        
    IL23_model_saved_chi2s.edges{end+1}=subconfiguration;
    IL23_model_saved_chi2s.chi2(end+1)=ar.chi2;
    IL23_model_saved_chi2s.p(end+1,:)=ar.p';
    
    end
    
    save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')
end

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

save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')
