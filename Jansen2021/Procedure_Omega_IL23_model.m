function [IL23_model_saved_chi2s]=Procedure_Omega_IL23_model(n_o,IL23_model_saved_chi2s)

edge_labels=IL23_model_saved_chi2s.edge_labels;

if n_o==length(edge_labels)
    disp('All significant edges are already part of the model.')
    return
end

%Retrieve current minimal model configurations:
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

disp(['Current found minimal model of size ',num2str(n_o),' has chi2=', num2str(IL23_model_saved_chi2s.chi2s(n_o)),', with edges:',])
IL23_model_saved_chi2s.initial_model{n_o}

added_edges=IL23_model_saved_chi2s.initial_model{n_o};
removed_edges=edge_labels(~contains(edge_labels,added_edges));

free_edges=[];
fixed_edges=[];
for m=1:n_o
    % We remove the edges of the n_o edge model one by one from the
    % N-sized model. If the removal of a single edge already results in a chi2
    % that is larger than the found chi2 of the n_o edge model, we know
    % this edge is essential for the n_o edge model and we fix it.
    disp('Compute subconfiguration:')
    subconfiguration=edge_labels(~contains(edge_labels,added_edges(m)))
    [chi2,IL23_model_saved_chi2s]=compute_chi2(subconfiguration,IL23_model_saved_chi2s);
    if chi2<IL23_model_saved_chi2s.chi2s(n_o)
        free_edges=[free_edges, added_edges(m)];
    else
        fixed_edges=[fixed_edges, added_edges(m)]; % We can't remove this edge
    end
end
if ~isempty(free_edges)
    % We enumerate all possible ways to remove combinations of the free
    % edges:
    MM={};
    for j=1:length(free_edges)
        MM{j}=nchoosek(free_edges,j);
    end
    
    for j=1:min(length(free_edges),length(removed_edges))
        M=MM{j};
        if length(removed_edges)-j>0
            P=nchoosek(removed_edges,length(removed_edges)-j);
        else
            P=[];
        end
        for m=1:size(M,1)
            % We check whether the removal of this combination of edges
            % from the N-sized model
            % results in a chi2 that is larger than the found chi2 of
            % the initial n_o-sized model. If so, this combination is essential and can not be removed and
            % we remove it from the list of combinations to try
            subconfiguration=edge_labels(~contains(edge_labels,M(m,:)));
            [chi2,IL23_model_saved_chi2s]=compute_chi2(subconfiguration,IL23_model_saved_chi2s);
            if chi2>IL23_model_saved_chi2s.chi2s(n_o)
                for q=j+1:length(free_edges)
                    temp=MM{q};
                    temp=temp(sum(contains(temp,M(m,:)),2)~=j,:);
                    MM{q}=temp;
                end
            else
                % We found a combination that might be replaced by edges
                % not present in the initial n_o edge model. We check
                % this:
                disp('We check whether the combination of edges')
                M(m,:)
                disp('could be replaced by other edges not part of the initial configuration.')
                
                for p=1:size(P,1)
                    disp('Check subconfiguration:')
                    subconfiguration=edge_labels(~contains(edge_labels,[M(m,:),P(p,:)]))
                    %Check whether this subconfiguration is contained
                    %in an already computed larger configuration:
                    flag=0;
                    for y=1:length(IL23_model_saved_chi2s.edges)
                        if sum(contains(IL23_model_saved_chi2s.edges{y},subconfiguration))==size(subconfiguration,2)
                            if IL23_model_saved_chi2s.chi2(y)>IL23_model_saved_chi2s.chi2s(n_o)
                                flag=1;
                            end
                        end
                    end
                    if flag==0
                        [chi2,IL23_model_saved_chi2s]=compute_chi2(subconfiguration,IL23_model_saved_chi2s);
                        if chi2<IL23_model_saved_chi2s.chi2s(n_o)
                            disp(['Found a better model configuration for model size ',num2str(n_o)])
                            %We repeat the procedure for the new model
                            %configuration:
                            Procedure_Omega_IL23_model(n_o,IL23_model_saved_chi2s)
                            return
                        end
                    end
                end
            end
        end
    end
end

%Set current minimal model configurations for every model size n:
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

%Set current relative AIC values for every model size n:
for i=1:size(IL23_model_saved_chi2s.initial_model,2)
    number_of_edges=size(IL23_model_saved_chi2s.initial_model{i},2);
    IL23_model_saved_chi2s.AIC(i)=IL23_model_saved_chi2s.chi2s(i)+2*number_of_edges;
end
save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')

end

function [chi2,IL23_model_saved_chi2s]=compute_chi2(subconfiguration,IL23_model_saved_chi2s)
edge_labels=IL23_model_saved_chi2s.edge_labels;

%Test if already computed:
flag=0;
for y=1:length(IL23_model_saved_chi2s.edges)
    if isequal(sort(unique(subconfiguration)),sort(unique(IL23_model_saved_chi2s.edges{y})))
        %disp('Already computed')
        chi2=IL23_model_saved_chi2s.chi2(y);
        flag=1;
    end
end
if ~flag
    arInit;
    ar=IL23_model_saved_chi2s.ar;
    
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
        arSetPars(subconfiguration{z},-2,1,[],-2,1);
    end
    
    arFitLHS_Jansen2021(51);
    
    if sum([((ar.p(ar.qFit==1)-ar.lb(ar.qFit==1))<=0.01),((ar.ub(ar.qFit==1)-ar.p(ar.qFit==1))<=0.01)])
        %disp('Warning: value of p outside of initial-guess sampling interval.')
    end
    
    if ar.LhsSampleSizeCalculation.D > 1
        IL23_model_saved_chi2s.multiple_minima_number_of_minima(end+1)=ar.LhsSampleSizeCalculation.D;
        IL23_model_saved_chi2s.multiple_minima_edges{end+1}=subconfiguration;
        IL23_model_saved_chi2s.multiple_minima_sample_size_calculation{end+1}=ar.LhsSampleSizeCalculation;
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
    chi2=ar.chi2;
    save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')
end
end
