%% Setup:
try
    load('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')
catch
    IL23_model_saved_chi2s.edge_labels=ar.pLabel(7:25);
    IL23_model_saved_chi2s.ar=ar;
end
edge_labels=IL23_model_saved_chi2s.edge_labels;
ar=IL23_model_saved_chi2s.ar;

%% No edges:
disp('Model with no edges:')
try  IL23_model_saved_chi2s.chi2(1);
    disp('Using Chi2 from IL23_model_saved_chi2s.mat.')
catch
    disp('Compute chi2:')
    arSetPars('LPS_IL10_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL1A_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL1B_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL23_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL6_c',0,1,[],-1.5,0);
    arSetPars('LPS_TNFA_c',0,1,[],-1.5,0);
    for i=1:length(edge_labels)
        arSetPars(edge_labels{i},-100,2,[],-100,1);
    end
    
    arFitLHS_Jansen2021(51);

    IL23_model_saved_chi2s.edges{1}={};
    IL23_model_saved_chi2s.chi2(1)=ar.chi2;
    IL23_model_saved_chi2s.p(1,:)=ar.p;
    IL23_model_saved_chi2s.multiple_minima_number_of_minima(1)=0;
    IL23_model_saved_chi2s.multiple_minima_number_of_minima_observed_once(1)=0;
    IL23_model_saved_chi2s.multiple_minima_edges{1}={};
   
    save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')
end

%% Full model:

disp('Compute chi2 for configurations with all significant edges):')
flag=0;
for y=1:length(IL23_model_saved_chi2s.edges)
    if max(size(IL23_model_saved_chi2s.edges{y}))==max(size(edge_labels))
        disp('Using Chi2 from IL23_model_saved_chi2s.mat.')
        flag=1;
    end
end

if ~flag
    
    arSetPars('LPS_IL10_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL1A_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL1B_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL23_c',0,1,[],-1.5,0);
    arSetPars('LPS_IL6_c',0,1,[],-1.5,0);
    arSetPars('LPS_TNFA_c',0,1,[],-1.5,0);
    for i=1:max(size(edge_labels))
        arSetPars(edge_labels(i),-1.5,1,[],-2,1); %Adds all edges
    end
    arFitLHS_Jansen2021(51)
    
    if sum([((ar.p(ar.qFit==1)-ar.lb(ar.qFit==1))<=0.01),((ar.ub(ar.qFit==1)-ar.p(ar.qFit==1))<=0.01)])
        %disp('Warning: returned a parameter value outside of initial-guess sampling interval.')
    end
    
    if ar.LhsSampleSizeCalculation.D > 1
        IL23_model_saved_chi2s.multiple_minima_number_of_minima(end+1)=ar.LhsSampleSizeCalculation.D;
        IL23_model_saved_chi2s.multiple_minima_edges{end+1}=edge_labels;
        IL23_model_saved_chi2s.multiple_minima_number_of_minima_observed_once(end+1)=ar.LhsSampleSizeCalculation.f1;
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
    
    IL23_model_saved_chi2s.edges{end+1}=edge_labels;
    IL23_model_saved_chi2s.chi2(end+1)=ar.chi2;
    IL23_model_saved_chi2s.p(end+1,:)=ar.p';
    
    save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')
end
