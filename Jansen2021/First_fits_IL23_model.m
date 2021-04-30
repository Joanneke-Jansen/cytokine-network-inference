%% Setup:
N=51; %Number of LHS fits
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
    
    arFitLHS_Jansen2021(N);

    IL23_model_saved_chi2s.edges{1}={};
    IL23_model_saved_chi2s.chi2(1)=ar.chi2;
    IL23_model_saved_chi2s.p(1,:)=ar.p;
    save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')
end

%% One edge:
for i=1:max(size(edge_labels))
    disp('Compute chi2 for model with edge:')
    disp(edge_labels{i})
    flag=0;
    for y=1:length(IL23_model_saved_chi2s.edges)
        if max(size(IL23_model_saved_chi2s.edges{y}))==max(size(edge_labels(i)))
            if sum(contains(IL23_model_saved_chi2s.edges{y},edge_labels{i}))==max(size(edge_labels(i)))
                disp('Using Chi2 from IL23_model_saved_chi2s.mat.')
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
        for l=1:length(edge_labels)
            arSetPars(edge_labels{l},-100,2,[],-100,1); %Removes all edges
        end
        arSetPars(edge_labels(i),-1.5,1,[],-2,1); %Inserts edge of interest
        arFitLHS_Jansen2021(51);
        
        if ar.LhsSampleSizeCalculation.D > 1
            if max(min(ar.ps(ar.chi2s-min(ar.chi2s)>1e-4,ar.qFit==1)'))>-16
                disp('More than one minimum found with all parameter values larger than 10^(-16).')
            return
            end
        end
        IL23_model_saved_chi2s.edges{end+1}=edge_labels(i);
        IL23_model_saved_chi2s.chi2(end+1)=ar.chi2;
        IL23_model_saved_chi2s.p(end+1,:)=ar.p';
        save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')
    end
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
    arFitLHS_Jansen2021(N)
        if ar.LhsSampleSizeCalculation.D > 1
            if max(min(ar.ps(ar.chi2s-min(ar.chi2s)>1e-4,ar.qFit==1)'))>-16
                disp('More than one minimum found with all parameter values larger than 10^(-16).')
            return
            end
        end
    IL23_model_saved_chi2s.edges{end+1}=edge_labels;
    IL23_model_saved_chi2s.chi2(end+1)=ar.chi2;
    IL23_model_saved_chi2s.p(end+1,:)=ar.p';
    
    save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')
end
