function [U, W, Pu, coeffMat, Pw] = samplePopulationISI(N, K, m, allPermutations, type, calculateWDistribution, calcPuAndNotPw)
    if ~exist('type', 'var')
        type = 'symmetric';
    end
    if ~exist('calculateWDistribution', 'var')
        calculateWDistribution = 1;
    end
    if ~exist('calcPuAndNotPw','var')
        calcPuAndNotPw = 0;
    end
    if m < K 
        fprintf('Notice! m < K, m should be greater than K!\n')
    end
    % ====== Generate infection coefficients matrix ======
    % P(Ui) = sum{ a_ij * P(Uj)}
    % a_ij ≠ 0 for m different items 
    
    if strcmp(type, 'symmetric')
    % Generate a (N div m) blocks matrix
    % each item is affected by its m-members family
    
        coeffMat = zeros(N,N);
        % fill N div m blocks
        curIdx = 1;
        while curIdx <= N-m
            for ii=1:m+1
                tempCoeffMat = sort(randfixedsum(m,1,1,0,1));
                relevantIdx = setdiff(curIdx:curIdx+m,[curIdx+ii-1]);
                coeffMat(curIdx + ii-1, relevantIdx) = tempCoeffMat;
            end
            curIdx = curIdx + m + 1;
        end
        % fill the rest mod(N, m) lines 
        while curIdx <= N
            tempCoeffMat = randfixedsum(m,1,1,0,1);
            relevantIdx = setdiff(curIdx:N, [curIdx]);
            anotherLegalAffected = 1:curIdx-1;
            relevantIdx = [relevantIdx, randsample(anotherLegalAffected , m-length(relevantIdx))];
            coeffMat(curIdx, relevantIdx) = tempCoeffMat;
            curIdx = curIdx + 1;
            
        end
    elseif strcmp(type, 'asymmetric')
    % each item is affected by the m-items before
        coeffMat = zeros(N,N);
        for ii=m+1:N
            coeffMat(ii, ii-m:ii-1) = sort(randfixedsum(m,1,1,0,1));
        end
        % fill the m-first rows in cyclic way
        for ii=1:m
            tempCoeffMat = randfixedsum(m,1,1,0,1);
            if length(mod(N+ii-m, N+1):mod(N+ii-1,N+1)) == m
                coeffMat(ii, mod(N+ii-m, N+1):mod(N+ii-1,N+1)) = tempCoeffMat(ii);
            else
                lastIdxList = mod(N+ii-m, N+1):N;
                coeffMat(ii, mod(N+ii-m, N+1):N) = tempCoeffMat(1:length(lastIdxList));
                numOfFirstIdxList = m - length(lastIdxList);
                coeffMat(ii, 1:numOfFirstIdxList) = tempCoeffMat(length(lastIdxList)+1:end);
            end
        end
    elseif strcmp(type, 'general')
        coeffMat = zeros(N,N);
        for ii=1:N
            coeffMat(ii, :) = randfixedsum(N,1,1,0,1); % random vector with length N and sum 1
        end
    end
    
    % ======= Sample population =======
    % choose the first defective item
    firstDefectiveItem = randi(N);
    [probabilityToDefective, U] = spreadInfectionUsingCorrMatrix(N, K, coeffMat, firstDefectiveItem);
    defectiveItems = find(U == 1);
    if calcPuAndNotPw
        W = -1; 
    else
        [tf,W] = ismember(defectiveItems, allPermutations, 'rows');
        if ~tf
            error('Error: permutation does not exist');
        end
    end
        
    % ======= Calculate W and U distribution if needed =======
    if calculateWDistribution && ~isempty(allPermutations)
        
        [Pw, Pu] = calculateIndicesStatistics(N, K, coeffMat, allPermutations);
        
    else
        if calcPuAndNotPw
            Pu = calculateUIndicesStatistics(N, K, coeffMat);
            Pw = [];
        else
            numOfW = size(allPermutations,1);    
            Pw = zeros(numOfW ,1);
            Pu = ones(N,1)/N;
        end
    end
    
if 0
    %% show that each item has the same probability to be detfective
    clear
    clc
    nmc = 1;
    N = 100;
    K = 10;
    m = 20; % if mod( N, m) ~= 0 you wont see the uniform distribution
    totU = zeros(1,N);
    Pw_tot = zeros(1,N);
    sumPu = zeros(1,N);
    errorCount = 0; % count error of W is not the index for U in the permutations array
    allPermutations = [];%getAllPossiblePermutations(N,K);
    calcPuAndNotPw = 1;
    for ii=1:nmc
        [U, W, Pu, coeffMat, Pw] = samplePopulationISI(N, K, m, allPermutations, 'asymmetric', 1, calcPuAndNotPw);
        totU = totU + U;
        sumPu(ii) = sum(Pu);
        if ~isempty(allPermutations) && ~isempty(setdiff(allPermutations(W,:),find(U==1)))
            errorCount = errorCount+1;
        end
    end
    if ~isempty(Pw)
        Pw_tot = Pw_tot + Pw;
    end
    totU = totU * 100 / nmc;
    figure;
    plot(1:N, totU(1:end));
    xlabel('item')
    ylabel('#time\_it\_was\_defective')
    ylim([0,100]);
    grid on;
    figure;imshow(coeffMat)
    colormap turbo
    
%     figure;histogram(sumPu)
    
end
end

function [probabilityToDefective, U] = spreadInfectionUsingCorrMatrix(N, K, coeffMat, firstDefectiveItem)
% sample population
% N - population size
% K - the number of the defective items
% coeffMat - NxN matrix where the i-th row includes the infection rate from
%           the i-th item to evryone else
% firstDefectiveItem - index of the first item that is known to be
%                       defective
% probabilityToDefective - Nx1 vector, the i-th element is the probaility
%                           of the i-th item to be defective after the 
%                           infection spread
% U - Nx1 boolean vector of the items. 1 = defective, 0 = not defective

    probabilityToDefective = zeros(N,1);
    probabilityToDefective(firstDefectiveItem,1) = 1;

    % spread the infection to the other items in steps
    currentInfectionStep = [firstDefectiveItem]; % a list of the defective 
                        % items that may affect others in this current step
    checkedItems = []; % list of items that could be affected
    while length(checkedItems) < N && ~isempty(currentInfectionStep)
        nextInfectionLevel = [];
        for ii=1:length(currentInfectionStep)
            currentDefectiveItem = currentInfectionStep(ii);
            % find all items the current item affect that were not checked
            % yet and are not going to be checked in the current infection
            % step
            relevantColumn = coeffMat(:,currentDefectiveItem);
            maybeInfectedList = find(relevantColumn ~= 0);
            maybeInfectedList = setdiff(maybeInfectedList, [checkedItems; currentInfectionStep]);
            probabilityToDefective(maybeInfectedList,1) = probabilityToDefective(maybeInfectedList,1) ...
                                                        + coeffMat(maybeInfectedList, currentDefectiveItem) ...
                                                        * probabilityToDefective(currentDefectiveItem,1);
            nextInfectionLevel = [nextInfectionLevel; maybeInfectedList];
        end 
        checkedItems = unique([checkedItems; currentInfectionStep]);
        currentInfectionStep = unique([nextInfectionLevel]);
    end
    
    % Choose the K items with the max probability
    [~, maxKIdx] = maxk(probabilityToDefective, K);
    U = logical(zeros(1, N));
    U(1, maxKIdx) = 1;
    
end

function [Pw, Pu] = calculateIndicesStatistics(N, K, coeffMat, allPermutations)
% N - population size
% K - the number of the defective items
% coeffMat - NxN matrix where the i-th row includes the infection rate from
%           the i-th item to evryone else
% allPermutations - nchoosek x K matrix in which each row is a possible
%                   permutation of K items in population of size N
%                   permutation is K indices of the defective items.
    numOfW = size(allPermutations,1);    
    Pw = zeros(numOfW ,1);
    Pu = zeros(N, 1);
    for ii=1:N
        % for each item as the initial defective item
        % spread the infection and get W
        firstDefectiveItem = ii;
        [probabilityToDefective_ii, U_ii] = spreadInfectionUsingCorrMatrix(N, K, coeffMat, firstDefectiveItem);
        defective_idx = find(U_ii == 1);
        Pu = Pu + U_ii';
        tf = ismember(allPermutations,defective_idx,'rows');
        w_idx = find(tf == 1);
        Pw(w_idx, :) = Pw(w_idx, :) + 1;
    end % TODO: verify that I should take the K deffecive items and not the mean probabilities to be defective 
    Pw = Pw / N;
    Pu = Pu / (N * K); % normalization to sum(Pu)=1
end

function  Pu = calculateUIndicesStatistics(N, K, coeffMat)
% N - population size
% K - the number of the defective items
% coeffMat - NxN matrix where the i-th row includes the infection rate from
%           the i-th item to evryone else
% allPermutations - nchoosek x K matrix in which each row is a possible
%                   permutation of K items in population of size N
%                   permutation is K indices of the defective items.
    Pu = zeros(N, 1);
    for ii=1:N
        % for each item as the initial defective item
        % spread the infection and get W
        firstDefectiveItem = ii;
        [probabilityToDefective_ii, U_ii] = spreadInfectionUsingCorrMatrix(N, K, coeffMat, firstDefectiveItem);
        defective_idx = find(U_ii == 1);
        Pu = Pu + U_ii';
    end % TODO: verify that I should take the K defecive items and not the mean probabilities to be defective 
    Pu = Pu / (N * K); % normalization to sum(Pu)=1
end
