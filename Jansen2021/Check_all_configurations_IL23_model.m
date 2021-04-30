function [IL23_model_saved_chi2s]=Check_all_configurations_IL23_model(I,IL23_model_saved_chi2s)

edge_labels=IL23_model_saved_chi2s.edge_labels;

%Compute chi2 for all models with one edge removed:
index=0;
min_chi2=IL23_model_saved_chi2s.chi2(1);
for w=1:length(edge_labels)
    subconfiguration=edge_labels(~contains(edge_labels,edge_labels(w)));
    [chi2, IL23_model_saved_chi2s]=compute_chi2(subconfiguration, IL23_model_saved_chi2s);
    if chi2<min_chi2
        index=w;
        min_chi2=chi2;
    end
end

IL23_model_saved_chi2s.optimal_model{length(edge_labels)-1}=edge_labels(~contains(edge_labels,edge_labels(index)));
IL23_model_saved_chi2s.chi2s(length(edge_labels)-1)=min_chi2;
save('IL23_model_saved_chi2s.mat','IL23_model_saved_chi2s')

%% Find optimal model configuration for size I:
if I>2 && I<length(edge_labels)-1
    disp(['Current found minimal model of size ',num2str(I),' has chi2=', num2str(IL23_model_saved_chi2s.chi2s(I)),', with edges:',])
    IL23_model_saved_chi2s.optimal_model{I}
    
    added_edges=IL23_model_saved_chi2s.optimal_model{I};
    removed_edges=edge_labels(~contains(edge_labels,added_edges));
    
    free_edges=[];
    fixed_edges=[];
    for m=1:I
        % We remove the edges of the N edge model one by one from the full
        % model. If the removal of a single edge already results in a chi2
        % that is larger than the found chi2 of the N edge model, we know
        % this edge is essential for the N edge model and we fix it.
        subconfiguration=edge_labels(~contains(edge_labels,added_edges(m)));
        [chi2,IL23_model_saved_chi2s]=compute_chi2(subconfiguration,IL23_model_saved_chi2s)
        if chi2<IL23_model_saved_chi2s.chi2s(I)
            free_edges=[free_edges, added_edges(m)]
        else
            fixed_edges=[fixed_edges, added_edges(m)] % We can't remove this edge
        end
    end
    if ~isempty(free_edges)
        % We enumerate all possible ways to remove combinations of the free
        % edges:
        MM={};
        for j=1:length(free_edges)
            MM{j}=nchoosek(free_edges,j)
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
                % from the full model
                % results in a chi2 that is larger than the found chi2 of
                % the initial core
                % N edge model. If so, this combination can not be removed and
                % we remove it from the list of combinations to try
                subconfiguration=edge_labels(~contains(edge_labels,M(m,:)));
                [chi2,IL23_model_saved_chi2s]=compute_chi2(subconfiguration,IL23_model_saved_chi2s);
                if chi2>IL23_model_saved_chi2s.chi2s(I)
                    for q=j+1:length(free_edges)
                        temp=MM{q};
                        temp=temp(sum(contains(temp,M(m,:)),2)~=j,:);
                        MM{q}=temp;
                    end
                else
                    % We found a combination that might be replaced by edges
                    % not present in the initial core N edge model. We check
                    % this:
                    disp('We check whether the combination of edges')
                    M(m,:)
                    disp('could be replaced by other edges not part of the initial configuration.')
                    
                    min_chi2=num2str(IL23_model_saved_chi2s.chi2s(I)); %The current minimal chi2;
                    chi2=0;
                    f=I;
                    try
                        while chi2<min_chi2
                            subconfiguration=[added_edges(~contains(added_edges,M(m,:))),edge_labels(~contains(edge_labels,IL23_model_saved_chi2s.optimal_model{f}))];
                            [chi2,IL23_model_saved_chi2s]=compute_chi2(subconfiguration,IL23_model_saved_chi2s);
                            f=f+1;
                        end
                    catch
                    end
                    
                    for p=1:size(P,1)
                        disp('Check subconfiguration:')
                        subconfiguration=edge_labels(~contains(edge_labels,[M(m,:),P(p,:)]))
                        %Check whether this subconfiguration is contained
                        %in an already computed larger configuration:
                        flag=0;
                        for y=1:length(IL23_model_saved_chi2s.edges)
                            if sum(contains(IL23_model_saved_chi2s.edges{y},subconfiguration))==size(subconfiguration,2)
                                if IL23_model_saved_chi2s.chi2(y)>IL23_model_saved_chi2s.chi2s(I)+1
                                    flag=1;
                                end
                            end
                        end
                        if flag==0
                            [chi2,IL23_model_saved_chi2s]=compute_chi2(subconfiguration,IL23_model_saved_chi2s);
                            if chi2<IL23_model_saved_chi2s.chi2s(I)
                                error(['Found a better model configuration for model size ',num2str(I)])
                                return
                            end
                        end
                    end
                end
            end
        end
    end
end
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
                error('More than one minimum found with all parameter values larger than 10^(-16).')
                return
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
