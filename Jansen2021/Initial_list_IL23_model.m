%% Starting from one edge, we add edges to the model to create an initial list of network configurations
disp('Starting from one edge, we add edges to the model to create an initial list of network configurations')
for n=2:length(edge_labels)
    
    %Of all already computed network configurations, find the best model of size n-1:
    index=0;
    min_chi_2=max(IL23_model_saved_chi2s.chi2);
    for y=1:length(IL23_model_saved_chi2s.edges)
        if size(IL23_model_saved_chi2s.edges{y},2)==(n-1)
            chi2=IL23_model_saved_chi2s.chi2(y);
            if chi2<min_chi_2
                index=y;
                min_chi_2=chi2;
            end
        end
    end
    IL23_model_saved_chi2s.selected{n-1}=IL23_model_saved_chi2s.edges{index};
    IL23_model_saved_chi2s.chi2s(n-1)=min_chi_2;
    
    added_edges=IL23_model_saved_chi2s.selected{n-1};
    remaining_edges=edge_labels(~contains(edge_labels,added_edges));
    
    for w=1:length(remaining_edges)
        subconfiguration=[added_edges, remaining_edges(w)];
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
            %Remove all edges:
            for p=1:length(edge_labels)
                arSetPars(edge_labels{p},-100,2,[],-100,[]);
            end
            %Insert edges of interest:
            for z=1:size(subconfiguration,2)
                arSetPars(subconfiguration(z),-1,1,[],[],[]);
            end
            %Set bounds on initial guess intervals:
            ar.lb(ar.qFit==1)=IL23_model_saved_chi2s.lb(ar.qFit==1);
            ar.ub(ar.qFit==1)=IL23_model_saved_chi2s.ub(ar.qFit==1);
            
            arFitLHS_Jansen2021(51);
            
            if sum([((ar.p(ar.qFit==1)-ar.lb(ar.qFit==1))<=0.01),((ar.ub(ar.qFit==1)-ar.p(ar.qFit==1))<=0.01)])
                arPrint
                disp('Warning: value of p outside of initial-guess sampling interval.')
                
                if ar.LhsSampleSizeCalculation.D > 1
                    if max(min(ar.ps(ar.chi2s-min(ar.chi2s)>1e-4,ar.qFit==1)'))>-16
                        disp('More than one minimum found with all parameter values larger than 10^(-16).')
                        return
                    end
                end
                
            end
            
            IL23_model_saved_chi2s.edges{end+1}=subconfiguration;
            IL23_model_saved_chi2s.chi2(end+1)=ar.chi2;
            IL23_model_saved_chi2s.p(end+1,:)=ar.p';
            
            save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')
        end
    end
end
