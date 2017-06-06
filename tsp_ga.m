
function varargout = tsp_ga(varargin)

% Initialize default configuration (do not change this part)
defaultConfig.xy          = 10*rand(50,2); % xy is the coordinate matrix
defaultConfig.dmat        = [];            % dmat is the distance matrix
defaultConfig.popSize     = 200;
defaultConfig.numGen      = 2000;
defaultConfig.crossProb   = 0.25;
defaultConfig.mutProb     = 0.5;
defaultConfig.eliteFract  = 0.02;

% Interpret user configuration inputs (do not change this part)
if ~nargin
    userConfig = struct();
elseif isstruct(varargin{1})
    userConfig = varargin{1};
else
    try
        userConfig = struct(varargin{:});
    catch
        error('Expected inputs are either a structure or parameter/value pairs');
    end
end

% Override default configuration with user inputs (do not change this part)
configStruct = get_config(defaultConfig,userConfig);

% Extract configuration
xy          = configStruct.xy;   % xy is the coordinate matrix
dmat        = configStruct.dmat; % dmat is the distance matrix
popSize     = configStruct.popSize;
numGen      = configStruct.numGen;
crossProb   = defaultConfig.crossProb;
mutProb     = defaultConfig.mutProb;
eliteFract  = defaultConfig.eliteFract;

if isempty(dmat)
    nPoints = size(xy,1);
    a = meshgrid(1:nPoints);
    dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),nPoints,nPoints);
end

% Verify Inputs (do not change this part)
[N,dims] = size(xy);
[nr,nc] = size(dmat);
if N ~= nr || N ~= nc
    error('Invalid XY or DMAT inputs!')
end
n = N; % make sure you do not use this varaible n for other puposes (e.g. for loop iteration)




%%%%%%%%%%%%%%%%% Initialize the Population %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% You don't need to change this part%%%%%%%%%%%%%%%%%%%%%%%
pop = zeros(popSize,n); % Matrix pop maintains the current population
pop(1,:) = (1:n);
for k = 2:popSize
    pop(k,:) = randperm(n);
end
%%%%%%%%%%%%%%%%% End of Population Initialization %%%%%%%%%%%%%%%%%%%%%%%%



totalDist = zeros(1,popSize); 
% totalDist is the vector of distances.Each element of this vector corresponds
% to the total distance (i.e. length) of a member of the population.



%% Starting GA iterations. In each iteration, a new generation is created %%%%%%
for iter = 1:numGen
    
    % Function calcToursDistances evaluates Each population member and 
    % calculates the total distance of each member
    totalDist = calcToursDistances(pop, popSize, dmat, n);
    
    
    
    % Elite selection: you should use the information in matrix totalDist 
    % to select a fraction eliteFract of the best members of the current
    % population pop. Keep these elite members in matrix elitePop.
    % Your elite selection code goes here:
    % ...
    % ...
    % ...
    numberOfElites = round(eliteFract*popSize);
    elitePop = zeros(numberOfElites, n);
    [sortedDist, indices] = sort(totalDist);
    for i = 1:numberOfElites
        row = pop(indices(i),:);
        elitePop(i,:) = row;
    end
    %%%%%%% end of elite selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Selection of new population: use totalDist to clacluate the fitness
    % and the cummulative probability distribution of the current population
    % pop. Then apply a roulette wheel selection to create a new population.
    % Keep this new population in matrix newPop.
    % Your roulette wheel selection code goes here:
    % ...
    % ...
    % ...
    fitnessTotal = 0;
    for i = 1:length(totalDist)
        fitnessTotal = fitnessTotal + (1/totalDist(i));
    end
    
    p = [];
    q = [];
    for i = 1:length(totalDist)
      p_i = (1/totalDist(i))/fitnessTotal;
      p = [p, p_i];
    end
    
    for i = 1:length(p)
        q_i = 0;
        for j = 1:i
            q_i = q_i + p(i);
        end
        q = [q, q_i];
    end
    
    listOfChosenChromosomes = [];
    %popSize = 48;
    for z = 1:popSize
        r = rand; %should be random
        if(r < q(1))
            %välj chromosome 1
            listOfChosenChromosomes = [listOfChosenChromosomes, 1];
            continue;
        end
        for i = 2:q
            if(r > q(i-1))
                if(r < q(i))
                    %välj chromosome i
                    listOfChosenChromosomes = [listOfChosenChromosomes, i];
                    break;
                end
            end
        end
    end
    newPop = [];
    for i = 1:listOfChosenChromosomes
        newPop = [newPop, totalDist(i)]; %maybe change this????
    end
    %%%%%%% end of roulette wheel selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % Update distance vector totalDist according to the selected population
    % newPop. Your code for updating totalDist goes here:
    totalDist = newPop;
    %%%%%% end of totalDist update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % Use the updated totalDist to calculate the new fitness values and 
    % cummulative probability distribution. Your code goes here:
    fitnessTotal = 0;
    for i = 1:length(totalDist)
        fitnessTotal = fitnessTotal + (1/totalDist(i));
    end
    p = [];
    q = [];
    for i = 1:length(totalDist)
      p_i = (1/totalDist(i))/fitnessTotal;
      p = [p, p_i];
    end
    
    for i = 1:length(p)
        q_i = 0;
        for j = 1:i
            q_i = q_i + p(i);
        end
        q = [q, q_i];
    end
    %%%%%% end of fitness and probabilty distribution update %%%%%%%%%%%%%%
    
    
    
    
    % Cross-over operator: implement the cross-over procedure
    % described in the home assignment. Use the cross-over
    % probability crossProb to obtain the parents that should
    % be selected purely random from the population newPop.
    % Your code goes here:
    % your cross-over code
    numberOfCrossovers = crossProb*popSize; %fixa avrudning till jämnt tal
    for i = 1:(numberOfCrossovers/2)
        r1 = randi(popSize);
        r2 = randi(popSize);
        if(r1 == r2)
            r2 = randi(popSize); %cheating...but it will work!!!!!
        end
        P1 = pop(r1,:); %TODO: change to newPop later
        P2 = pop(r2,:); %TODO: change to newPop later
        T1 = P1;
        T2 = P2;
        K = floor(0.3*length(P1));
        e = '';
        O = [];
        
        %While Tl ? ? do
        while(length(T1) > 0)
        %    Split Tl in two sub-tours T11 and T12 such that length of T11 = min {length of T1, K};
            T11_length = length(T1);
            if(T11_length > K)
                T11_length = K;
            end
            T11 =  T1(1:T11_length);
            T12 =  T1(T11_length+1:length(T1));
        %    Append T11 to O;
            O = [O T11];
        %    Update Tl by removing from it the cities in T11
            T1 = T12;
        %    Update T2 by removing from it the cities in T11
            %?????????????????????
            for i = 1:length(T11)
                k = find(T2==T11(i));
                for j = 1:length(k)
                   T2(k(j)) = []; 
                end
            end
        %    X = Tl;
            X = T1;
        %    Tl = T2;
            T1 = T2;
        %    T2 = X;
            T2 = X;
        %End {while};
        end
        %Return O as the offspring
    end
    %%%%%%% End of cross-over operator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%% Mutation Operator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = rand();
    if r <= mutProb
        off_indx = randi([1 popSize], 1, 1); 
        offspring = pop(off_indx, :); % replace this line of code so that the offspring
                                      % will be the one created by the
                                      % cross-over operation.
       
        routeInsertionPoints = sort(ceil(n*rand(1,2)));
        I = routeInsertionPoints(1);
        J = routeInsertionPoints(2);
        
        % 2-opt mutation (simply swaps two cities)
        offspring([I J]) = offspring([J I]);
        
        % now, you should replace one of the parents invloved in
        % cross-over with this mutated offspring, then update the
        % population newPop.
    end
    
    %%%%%%%%%% End of mutation operator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % Now, it is time to replace the worst members of newPop with the elite 
    % members you stored in matrix elitePop (see Elite selection in the begining
    % of this iteration).
    % Your code goes here:
    % ...
    % ...
    % ...
    % ...
    %%%%%%% End of elite replacement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Finally, the new population newPop should become the current population.
    % pop = newPop;    % Uncomment this line when you finished all previous
                       % steps.

end
%%%%%% End of GA ietartions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now, we find the best route in the last generation (you don't need to
% change this part). The best route is returned in optRoute, and its 
% distance in minDist.
totalDist = calcToursDistances(pop, popSize, dmat, n);
[minDist,index] = min(totalDist);
optRoute = pop(index,:);

% Return Output (you don't need to change this part)
if nargout
    resultStruct = struct( ...
        'optRoute',    optRoute, ...
        'minDist',     minDist);
    
    varargout = {resultStruct};
end

end


% The following code is for configuation of user input (do not change
% this). Subfunction to override the default configuration with user inputs
function config = get_config(defaultConfig,userConfig)

% Initialize the configuration structure as the default
config = defaultConfig;

% Extract the field names of the default configuration structure
defaultFields = fieldnames(defaultConfig);

% Extract the field names of the user configuration structure
userFields = fieldnames(userConfig);
nUserFields = length(userFields);

% Override any default configuration fields with user values
for i = 1:nUserFields
    userField = userFields{i};
    isField = strcmpi(defaultFields,userField);
    if nnz(isField) == 1
        thisField = defaultFields{isField};
        config.(thisField) = userConfig.(userField);
    end
end

end

