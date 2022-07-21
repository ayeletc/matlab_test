clear
%% Count #possiblyDefected after CoMa and DD
% 1.
% 1.1. CoMa with T=Tml
% 1.2. count PD1 (should be ~2k)
% 2.
% 2.1. DD 
% 2.2. count PD2

%% Config simulation
N               = 100;
vecK            = 2:4:20;%1:20:100;%1:5:30;%1:20:150; %10;%round(beta * N ^ alpha);
sampleMethod    = 'ISI'; 
isiType         = 'asymmetric';
m               = 20;
Pemax_theor     = 0.01;%0.it 05:0.1:0.95; % for N = 500, K = 5, possible T range for DND+MAP: (35, ∞)
nmc             = 500;
saveRaw         = 0; 
savePath        = '/Users/ayelet/Library/CloudStorage/OneDrive-Technion/Alejandro/count_possibly_defected_results/count_pd_and_dd/';
isPlot          = 1;
doMAP           = 0;
rng(123)
allPermutations = [];
calculateWDistribution = 1;
calcOnlyPu = 1;
enlargeTestsNumByFactors = [0.5 0.75 1 1.5];%[0.5 1 1.5];%[0.25 1 2];
Tbaseline = 'ML'; % options: 'ML', 'lb_no_priors', 'lb_with_priors'
method_DD = 'Normal'; % options: Normal, Iterative, Sum
invalid = -1;
%% Initialize counters
countSuccess = 0;
numOfK = length(vecK);
numOfTestScale = length(enlargeTestsNumByFactors);
count_DND1 = zeros(numOfK, numOfTestScale);
count_PD1 = zeros(numOfK, numOfTestScale);
count_DD2 = zeros(numOfK, numOfTestScale);
count_DND3 = zeros(numOfK, numOfTestScale);
count_PD3 = zeros(numOfK, numOfTestScale);
% vecT = ceil(1.01 * vecK * log(N)/log(2) .* enlargeTestsNumByFactor); % ceil((1-Pemax_theor) * K * log(N/K)); % try 
%% Start simulation

for idxK=1:numOfK
    K = vecK(idxK);
    fprintf('K = %d\n', K);
    % For each K calculate number of test according the Tml and scale factor
    if strcmp(Tbaseline, 'ML') 
        vecT = ceil(1.01 * K * log(N)/log(2) .* enlargeTestsNumByFactors); % ceil((1-Pemax_theor) * K * log(N/K));
    elseif ctrcmp(Tbaseline, 'lb_no_priors')
        vecT = ceil((1-Pemax_theor) * K * log(N/K) .* enlargeTestsNumByFactors);
    end
    for idxT=1:numOfTestScale
        T = vecT(idxT);
        for nn=1:nmc
    %         fprintf('nn = %d\n', nn);
            %% Sample
            [U, W, Pu, coeffMat, Pw] = samplePopulationISI(N, K, m, allPermutations, isiType, calculateWDistribution, calcOnlyPu);
            %% 1. Definitely Not Defective
            % Encoder - bernoulli 
            p = 1-2^(-1/K); % options: 1/K, log(2)/K, 1-2^(-1/K)
            X = rand(T, N) < p; % iid testing matrix
            testedMat = bsxfun(@and, X, U); 
            Y = sum(testedMat, 2) > 0;

            % Decoder - CoMa
            PD1 = 1:N;
            DND1 = [];
            for ii = 1:T
                if length(PD1) <= K
                    break 
                end
                if Y(ii) == 0
                    for jj = PD1
                        if X(ii,jj) == 1 % definitely not defected
                            PD1 = PD1(PD1 ~= jj);
                            DND1 = [DND1; jj];
                        end
                    end
                end
            end
            count_DND1(idxK, idxT) = count_DND1(idxK, idxT) + length(DND1);
            count_PD1(idxK, idxT) = count_PD1(idxK, idxT) + length(PD1);

            if length(PD1) <= K % all the PD are DD - all defective found
                count_DD2(idxK, idxT) = count_DD2(idxK, idxT) + length(PD1);
                continue
            end
            %% 2. WORKING Definite Defective
            % steps 1&2
            if strcmp(method_DD, 'Normal')
                DD2 = [];

                for ii = 1:T
                    if Y(ii) == 1 && sum(X(ii,PD1)) == 1 % only 1 item among the PD equals 1 and the rest equal 0
                        jj = find(X(ii,PD1) == 1); % find the definite defective item; index in PD1 array
                        defective = PD1(jj);
                        if ~sum(ismember(DD2, defective)) % add jj only if jj is not already detected as DD
                            DD2 = [DD2, defective ];
                        end
                    end
                end
                count_DD2(idxK, idxT) =  count_DD2(idxK, idxT) + length(DD2);

                if length(DD2) >= K % all defective found
                    fprintf('All defective found\n');
                    continue 
                end

            %% 2. TRY iterative Definite Defective - currently there is a bug:
            % TODO = check in the notes!
            % I need to see in my notes what we said about participating
            % item and sum=1
            % steps 1&2
            elseif strcmp(method_DD, 'Iterative')
                DD2 = [];
                tryAgainDD = 1;
                while tryAgainDD
                    tryAgainDD = 0;
                    for tt = 1:T
                        if Y(tt) == 0
                            continue
                        end
%                         X(tt,PD1)
                        % calculate sum over not DD
                        countParticipantWhoAreNotDD = 0;
                        participants = find(X(tt,PD1) == 1); % the indices of the PD that participate int the ii test
                                                            % the corresponding index in X = PD1(participants)
                        participantsWhoAreNotDD = ~ismember(PD1(participants), DD2);
                        countParticipantWhoAreNotDD = sum(participantsWhoAreNotDD);
                        if countParticipantWhoAreNotDD == 1 %&& participantsWhoAreNotDD % only 1 item among the PD (who have not already been detected as DD) equals 1 and the rest equal 0
                            jj = participants(participantsWhoAreNotDD); % find the definite defective item; index in PD1 array
                            defective = PD1(jj); % index of the defective in X (index between 1 to N)
                            if ~sum(ismember(DD2, defective)) % add jj only if jj is not already detected as DD
                                DD2 = [DD2, defective ];
                                tryAgainDD = 1;
                            end
                        end
                    end
                end
                count_DD2(idxK, idxT) =  count_DD2(idxK, idxT) + length(DD2);
                if length(DD2) >= K % all defective found
                    fprintf('All defective found\n');
                    continue 
                end
            %%
            elseif strcmp(method_DD, 'Sum')
                Y_sum = sum(testedMat, 2); % keep the levels
                iter = 0;
                helpful_iter = [];
                DD2 = [];
                DD_rows = [];
                tryAgainDD = 1;
                while tryAgainDD
                    tryAgainDD = 0;
                    iter = iter + 1;
                    for tt = 1:T
                        if Y(tt) == 0 || ismember(tt, DD_rows)% skip rows that based on them we detected DD
                            continue
                        end
                        countParticipantWhoAreNotDD = 0;
                        participants = find(X(tt,PD1) == 1); % the indices of the PD that participate int the ii test
                                                            % the corresponding index in X = PD1(participants)
                        participantsWhoAreNotDD = ~ismember(PD1(participants), DD2);
                        countParticipantWhoAreNotDD = sum(participantsWhoAreNotDD);
                        participantsWhoAreDD = ismember(PD1(participants), DD2);
                        countParticipantWhoAreDD = sum(participantsWhoAreDD);
                        if countParticipantWhoAreNotDD == 1 && Y_sum(tt)-countParticipantWhoAreDD == 1%&& participantsWhoAreNotDD % only 1 item among the PD (who have not already been detected as DD) equals 1 and the rest equal 0
                            jj = participants(~ismember(participants, participantsWhoAreDD));
                            defective = PD1(jj); % index of the defective in X (index between 1 to N)
                            if ~sum(ismember(DD2, defective)) % add jj only if jj is not already detected as DD
                                DD2 = [DD2, defective ];
                                DD_rows = [DD_rows, tt];
                                tryAgainDD = 1;
                                helpful_iter = [helpful_iter, iter];
                            end
                        end
                    end
                end
                count_DD2(idxK, idxT) =  count_DD2(idxK, idxT) + length(DD2);
                if length(DD2) >= K % all defective found
                    fprintf('All defective found\n');
                    continue 
                end
                if ~isempty(helpful_iter) && length(helpful_iter) > 1
                    fprintf([num2str(length(helpful_iter)) 'helpful iter\n']);
                end
            end
            
            % find all unknown
            unknown2 = PD1(~ismember(PD1, DD2));
            notDetectedDefectives = K-length(DD2);
            %% 3. MAP?
            % 3.1 calcualte new priors based on X(Y==1, PD1)
            if ~doMAP
                continue;
            end
            Pu3 = Pu;
            for tt=1:size(X,1)
                if Y(tt) == 0
                    continue
                end
                % calculate sum over not DD
                participants = find(X(tt,PD1) == 1); % the indices of the PD that participate int the ii test
                                                    % the corresponding index in X = PD1(participants)
                participantsWhoAreDD = ismember(PD1(participants), DD2);
                countParticipantWhoAreDD = sum(participantsWhoAreDD);
                if countParticipantWhoAreDD > 0
                    % what to do id there is a DD in the line?
                    continue;
                end
                % if there is no DD in the line - we give the participants priors
                nParticipants = length(participants);
                prior = 1+1/nParticipants; 
                participantsTrueIdx = PD1(participants);
                for participant=participantsTrueIdx
                    Pu3(participant) = Pu3(participant) * prior;
                end                
            end
            % normalize new priors:
            Pu3 = Pu3 / sum(Pu3);
                   
            % 3.2 start MAP
            allPermutations3 = nchoosek(1:length(unknown2),notDetectedDefectives); % not the correct indices
            allPermutations3 = unknown2(allPermutations3); % the correct indices in range[1,N]
            numOfPermutations3 = size(allPermutations3,1);
            
            apriori = invalid*ones(numOfPermutations3,1);
%             numOfPossibeComb = 0;
            for comb=1:numOfPermutations3
                % calculate Y for the w-th permutation 
                permute = allPermutations3(comb,:);
                U_forW = zeros(1,N);
                U_forW([permute DD2]) = 1;
                
                X_forW = bsxfun(@and, X, U_forW);
                Y_forW = sum(X_forW, 2) > 0;

                % possible case: Yw = Y
                % the case: Yw = 0 and Y = 1 is possible when T is too small
                % for example: U = [1 0 0]; T=1; X = [0 0 0]

                % skip the case: Yw = 1 and Y = 0:
                if ~isempty(find(Y_forW ~= Y))
                    continue
                end    
                % numOfPossibeComb = numOfPossibeComb + 1;
                apriori(comb) = 1;
                for tt=1:T
%                     if Y(tt) == 0 
%                         continue
%                     end
                    % calcualte P(Xsw(t))
                    % TODO = check in the notes!
                    numOfOnesInXwt = sum(X_forW(tt,:) == 1);
                    P_q = p^numOfOnesInXwt;
                    P_Xsw_t = P_q;
                    for ii=permute % TODO: do I need to add DD2?
                        P_Xsw_t = P_Xsw_t * Pu3(ii);
                    end
                    apriori(comb) = apriori(comb) * P_Xsw_t; %get zeros
                end
            end
            
            [maxAPriori, maxLikelihoodW] = max(apriori);
            estU = zeros(size(U));
            estU(allPermutations3(maxLikelihoodW,:)) = 1;  
            estU(DD2) = 1;
            if sum(U ~= estU) == 0
                countSuccess = countSuccess + 1;
            end

        end
    end
end
%% Normalize counters
count_DND1 = count_DND1 / nmc;
count_PD1 = count_PD1 / nmc;
count_DD2 = count_DD2 / nmc;
count_DND3 = count_DND3 / nmc;
count_PD3 = count_PD3 / nmc;

%% Visualize
markerStyles = ['o', 'p', 's','.'];
curveStyles = [':','--','!','-'];
curveColors = ['r','g','b','p','c','k'];
if isPlot
    figure;
    subplot(131)
    legendCell = {};
    for idxT=1:numOfTestScale
        
        scatter(vecK, count_DND1(:,idxT), 'filled', curveColors(idxT),markerStyles(idxT));
        hold on
        legendCell{end+1} = ['T=' num2str(vecT(idxT))];
    end
    hold off
    legend(legendCell);
    title({'#DND vs. K after CoMa', ['N = ' num2str(N), ', T=T_{Tbaseline}*[' num2str(enlargeTestsNumByFactors) ']'], ['#iterations=' num2str(nmc)]}, 'FontSize', 16)
    xlabel('K', 'FontSize', 16)
    ylabel('#DND(1)', 'FontSize', 16)
    ylim([0, N])
    grid on
    
    % Like Fig.2 in multi_level GT paper:
    subplot(132)
    legendCell = {};
    for idxT=1:numOfTestScale
        T = vecT(idxT);
        scatter(vecK, count_PD1(:,idxT), 'filled', curveColors(idxT), markerStyles(1));
        hold on
        theoreticalPD = vecK+(N-vecK).*((1-p*(1-p).^vecK).^T);   
        scatter(vecK, theoreticalPD, 'filled',curveColors(idxT), markerStyles(2));
        legendCell{end+1} = ['Empirical, T=' num2str(T)];
        legendCell{end+1} = ['Theoretical, T=' num2str(T)];
    end
    plot(vecK, vecK, [curveColors(end),curveStyles(3)])
    legendCell{end+1} = 'K';
    hold off
    legend(legendCell);
    title({'#PD vs. K after CoMa', ['N = ' num2str(N), ', T=T_{Tbaseline}*[' num2str(enlargeTestsNumByFactors) ']'], ['#iterations=' num2str(nmc)]}, 'FontSize', 16)
    xlabel('K', 'FontSize', 16)
    ylabel('#PD(1)', 'FontSize', 16)
    ylim([0, N])
    grid on
    
    % New:
    subplot(133)
    legendCell={};
    for idxT=1:numOfTestScale
        T = vecT(idxT);
        plot(vecK, count_DD2(:,idxT), [curveStyles(1), markerStyles(1), curveColors(idxT)], 'MarkerFaceColor',curveColors(idxT));
        hold on;
        plot(vecK, count_PD1(:,idxT) - count_DD2(:,idxT), [curveStyles(1), markerStyles(2), curveColors(idxT)], 'MarkerFaceColor',curveColors(idxT));
        legendCell{end+1} = ['DD(2), T=' num2str(T)];
        legendCell{end+1} = ['Unknown, T=' num2str(T)]; 
    end
    plot(vecK, vecK, [curveStyles(end), curveColors(end)])
    legendCell{end+1} = 'K';
    hold off
    legend(legendCell);
    title({'#DD vs. K after CoMa and DD', ['N = ' num2str(N), ', T=T_{' Tbaseline '}*[' num2str(enlargeTestsNumByFactors) ']'], ['#iterations=' num2str(nmc)]}, 'FontSize', 16)
    xlabel('K', 'FontSize', 16)
    ylabel('#DD(2)', 'FontSize', 16)
    ylim([0, N])
    grid on
    sgtitle([method_DD ' DD'], 'FontSize', 18)
end
%% Save
if saveRaw
    formatDateOut = 'ddmmyy_HHMM';
    timestr = datestr(now,formatDateOut);
%     vecTStr = strrep(strrep(num2str(orgVecT), ' ', '_'), '__','_');
    experimentStr = ['N' num2str(N) '_nmc' num2str(nmc)];
    save([savePath filesep 'countPDandDD_' experimentStr '_' timestr])
end