% initialize start position
start = [1 1];
% initialize finish position
finish = [500 500];
% Read a binary map from an image
map = im2bw(imread('random_map.bmp'));
% User selects a option for selection
selection = option1();
% User selects a option for crossover
crossOver = option2();
% User selects a option for mutation
mutation = option3();
numberofPoints = 10;
% This creates a population size (how many chromosomes in a population)
population_size = 1000;
% This defines how iteration the population is going to go through the Genetic algorthim
iter = 1000;
% This creates a matrix for the population with the size of population_size by 20 (the number of x and y coords)
population = zeros(population_size,numberofPoints*2);
% This creates a matrix for all the fittest chromosomes in each iteration
fittest = zeros(iter, 1);
% Starts a timer so we can see how long it takes
tic;


% This for loop iterates i from 1 to population_size and in the for loop there is a second for loop which
% increments j from 1 to 10. In the second for loop it will generate a random x and y coord between 1 and 500
% whole numbers only. Then it will check if its position on the binary image is a zero (shape) if it is then it will 
% keep on generating a random x and y coords till it is not zero. Then it store the x coord in population in position 
% i, j*2 - 1 and store the y coords in population in position i,j*2 it does this so x and y are always next to each other in
% in the matrix like (x,y,x,y,x,y).
for i = 1:population_size
    % Generate 10 random (X, Y) pairs within the map boundaries
    for j = 1:numberofPoints
        randX = randi(500,1);
        randY = randi(500,1);
        
        % Ensure the randomly selected position is within a valid region on the map
        while(map(randX,randY) == 0)
            randX = randi(500,1);
            randY = randi(500,1);
        end
        
        % Store the (X, Y) pairs in the population matrix
        population(i,j*2 - 1) = randX;
        population(i,j*2) = randY;
    end 
end
% Add a column of zeros to the population matrix for the fitness of each chromosome
population = [population zeros(population_size,1)];

% This for loop iterates k from 1 to iter then it enters anothe for loop which iterates i from 1 to population_size. 
% In the second for loop it will call the function Fitness_Function with the parameters population(i,:) and map.
% This function generates a fitness for the population a higher fitness the better and stores it in position i, size(population,2).
% Then it will do this for all the population_size. Then outside the second for loop it will sort the rows for population based on 
% its last column in ascending order. Then it will select the fittest and put it in position k,1. Then it will 
for k = 1:iter
    % calaculate the fitness for each chromosome in the population
    for i = 1:population_size
        population(i, size(population,2)) = Fitness_Function(population(i,:),map);
    end 
    % sort the population based on the fitness it ascending order
    population = sortrows(population,size(population,2));
    % pick the fittest chromosome and put it in the fittest matrix at position k,1
    fittest(k, 1) = population(end, size(population,2));
    % creates a new a population with the size of population_size by 20
    population_new = zeros(population_size,numberofPoints*2);
    % This selects the top 20% of the current population and copies them to the beginning of the new population
    population_new(1:(0.2*population_size),:) = population(population_size-(0.2*population_size-1):population_size,1:numberofPoints*2);
    % The sets the variable population_new_num to how many chromosomes where copied over from the current population
    % which is 20%
    population_new_num = (0.2*population_size);
    while (population_new_num < population_size)
        if selection == 0
            % calculates the weights for each member in the population by deviding their fitness by the 
            % sum of the population the fitness
            weights = population(:,(numberofPoints*2) + 1)/sum(population(:,(numberofPoints * 2) + 1));
            % selects a index by calling RouletteWheelSelection with the weights we just calaculated
            choice1 = RouletteWheelSelection(weights);
            % selects a index by calling RouletteWheelSelection with parameter weights
            choice2 = RouletteWheelSelection(weights);
        end
        if selection == 1
            % selects a index by calling TournamentSelection with the parameters 20 and population
            choice1 = TournamentSelection(20,population);
            % selects a index by calling TournamentSelection with the parameters 20 and population
            choice2 = TournamentSelection(20, population);
        end
        if selection == 2
            % creates a vector called ranks with integers from 1 to population_size
            ranks = 1:population_size;
            % weights is equal to ranks devided by the sum of the population fitness
            weights = ranks / sum(1:population_size);
            % selects a index by calling RouletteWheelSelection with the weights we just calaculated
            choice1 = RouletteWheelSelection(weights);
            % selects a index by calling RouletteWheelSelection with the weights we just calaculated
            choice2 = RouletteWheelSelection(weights);
        end
        % selects a chromosome based on the choice1
        temp_chromosome_1 = population(choice1, 1:(numberofPoints*2));
        % selects a chromosome based on the choice2
        temp_chromosome_2 = population(choice2, 1:(numberofPoints*2));
        % If crossover propability is less than 0.6 percent
        if (rand < 0.6)
            if crossOver == 0
            % calls the kpointCrossover function with the parameters temp_chromosome_1 and temp_chromosome_2 and 1
            % which returns the new temp_chromosome_1 and temp_chromosome_2 after they have been through kpointCrossover
                [temp_chromosome_1, temp_chromosome_2] = kpointCrossover(temp_chromosome_1, temp_chromosome_2, 1);
            end
            if crossOver == 1
            % calls the UniformCrossover function with the parameters temp_chromosome_1 and temp_chromosome_2 and 0.5
            % which returns the new temp_chromosome_1 and temp_chromosome_2 after they have been through UniformCrossover
                [temp_chromosome_1, temp_chromosome_2] = UniformCrossover(temp_chromosome_1, temp_chromosome_2, 0.5);
            end
        end
        % If the mutation propability is less than 0.4
        if (rand < 0.6)
            if mutation == 0
                % calls the Bit flip mutation mutation with the parameters temp_chromosome_1 and 20
                % which returns the temp_chromosome_1 after they have been
                % through bit flip
                temp_chromosome_1 = BitFlipMutation(temp_chromosome_1);
            end 
            if mutation == 1
                % calls the ValueEncodingMutation mutation with the parameters temp_chromosome_1 and 20
                % which returns the temp_chromosome_1 after they have been through ValueEncodingMutation
                temp_chromosome_1 = ValueEncodingMutation(temp_chromosome_1,50);
            end 
        end
        % If the mutation propability is less than 0.4
        if (rand < 0.6)
            if mutation == 0
                % calls the Bit flip mutation mutation with the parameters temp_chromosome_1 and 20
                % which returns the temp_chromosome_1 after they have been
                % through bit flip
                temp_chromosome_2 = BitFlipMutation(temp_chromosome_2);
            end 
            if mutation == 1
                % calls the ValueEncodingMutation mutation with the parameters temp_chromosome_2 and 20
                % which returns the temp_chromosome_2 after they have been through ValueEncodingMutation
                temp_chromosome_2 = ValueEncodingMutation(temp_chromosome_2,50);
            end 
        end
        % initializes tC1 and tC2 to 0
        tC1 = 0;
        tC2 = 0;
        for i = 1:numberofPoints*2
            % checks if temp_chromosome_1 at position i is greater than 500 or less than or equal to 0 (out of the map)
            if (temp_chromosome_1(i) > 500 || temp_chromosome_1(i) <= 0)
                % sets tC1 to 1
                tC1 = 1;
            end
            % checks if temp_chromosome_2 at position i is greater than 500 or less than or equal to 0 (out of the map)
            if (temp_chromosome_2(i) > 500 || temp_chromosome_2(i) <= 0)
                % sets tC2 to 1
                tC2 = 1;
            end
        end 
        if tC1 == 0
            % adds one to population_new_num
            population_new_num = population_new_num + 1;
            % adds temp_chromosome_1 to population_new at position population_new_num
            population_new(population_new_num,:) = temp_chromosome_1;
        end 
        if tC2 == 0
            % checks if population_new_num is smaller than population_size
            if (population_new_num < population_size)
                % adds one to population_new_num
                population_new_num = population_new_num + 1;
                % adds temp_chromosome_2 to population_new at position population_new_num
                population_new(population_new_num,:) = temp_chromosome_2;
            end
        end 
    end
    % sets population equal to population_new
    population(:,1:numberofPoints*2) = population_new;
end 
% stops the timer and sets its value to the variable timer
timer = toc;
% display how long the process took
disp("The process took " + timer)
% calculates the fitness for each member in the population
for i = 1:population_size
    population(i, size(population,2)) = Fitness_Function(population(i,:),map);
end 
% sort the population based on its fitness in ascending order
population = sortrows(population,(numberofPoints*2) + 1);
% get the solution by getting the last chromosome in the population
solution = (population(end, 1:numberofPoints*2));
% creates a figure and draws the solution on map
figure()
path = [start; [solution(2:2:end)' solution(1:2:end)']; finish];
clf;
imshow(map);
rectangle('position',[1 1 size(map)-1],'edgecolor','k');
line(path(:,2),path(:,1));
for i = 2:(numel(solution)/2) + 1
    text(path(i,2), path(i,1), num2str(i - 1), 'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
end

% This function requires a chromosome and the map
function fitness = Fitness_Function(tempChromosome, map)
    %tempChromosome
    % sets pathDistance to zero
    pathDistance = 0;
    % creates a chromosome with the size of 1 by 24
    chromosome = zeros(1,size(tempChromosome,2) + 4);
    % sets chromosome at position 1 to 1 (start x position)
    chromosome(1) = 1;
    % sets chromosome at position 2 to 1 (start y position)
    chromosome(2) = 1;
    % sets chromosome at position 23 to 500 (end x position)
    chromosome(size(tempChromosome,2) + 2) = 500;
    % sets chromosome at position 24 to 500 (end y position)
    chromosome(size(tempChromosome,2) + 3) = 500;
    % sets the values for chromosome at position 3 to 22 equal to the values at tempChromosome 
    % 1 to 20
    chromosome(3:size(tempChromosome,2) + 1) = tempChromosome(1:size(tempChromosome,2) - 1);
    % sets the intersection to 0
    intersection = 1;
    % this iterates i from 1 to size(chromosome,2)-3 by 2 each time
    for i = 1:2:size(chromosome,2)-3
        % sets check equal to 0
        check = 0;
        % sets checkNaN equal to 0
        checkNaN = 0;
        % x1 is equal to chromosome at position i
        x1 = chromosome(i);
        % x2 is equal to chromosome at position i + 2
        x2 = chromosome(i+2); 
        % y1 is equal to chromosome at position i + 1
        y1 = chromosome(i+1); 
        % y2 is equal to chromosome at position i + 3
        y2 = chromosome(i+3);
        % This adds the euclidian distance between the two points x1,y1 and x2,y2 to the variable pathDistance
        pathDistance = pathDistance + sqrt((x2-x1)^2 + (y2 - y1)^2);
        % This calaculates the gradient of the two points x1 y1 and x2,y2
        gradient = ((y2 - y1) / (x2 - x1));
        % If the gradient is equal to NaN then x1 = x2 = y1 = y2
        if isnan(gradient)
            % sets checkNaN to 1
            checkNaN = 1;
            % then checks if the point on the map is equal to 0 (shape)
            if(map(y1,x1) == 0)
                % adds a penalty to the chromosome
                intersection = intersection + 2000;
            end
        end
        % If the gradient is equal to inf or -inf and checkNaN == 0 then the line is a vertical line
        if gradient == inf || gradient == -inf & checkNaN == 0
            % sets check to 1
            check = 1;
            % sets opt to 0
            opt = 0;
            % calaculates the startPos and finishPos of the vertical line using calcOrder(y1,y2)
            [startPos, finishPos] = calcOrder(y1,y2);
            % this calculates the penalty of the vertical line using the function straightLine with the parameters 
            % startPos, finishPos, x1, opt, intersection, map
            intersection = straightLine(startPos, finishPos, x1, opt, intersection, map);
        end
        % If the gradient is equal to 0 and checkNaN is 0 then the line is a horizontal line
        if gradient == 0 & checkNaN == 0
            % sets check to 1
            check = 1;
            % sets opt to 1
            opt = 1;
            % calaculates the startPos and finishPos of the horizontal line using calcOrder(x1,x2)
            [startPos, finishPos] = calcOrder(x1,x2);
            % this calculates the penalty of the horizontal line using the function straightLine with the parameters 
            % startPos, finishPos, y1, opt, intersection, map
            intersection = straightLine(startPos, finishPos, y1, opt, intersection, map);
        end
        % this calaculates the y intercept
        c = (gradient * (-x1)) + y1;
        % if the absolute value of x1 + x2 is less than or equal to the absolute value of y1 and y2 and check is zero
        % and checkNaN is zero then this will use the go along the x axis to calaculate the y value
        if abs(x1 + x2) <= abs(y1 + y2) & check == 0 & checkNaN == 0
            % sets opt to 0
            opt = 0;
            % calculates the startPos and finishPos of the line using calcOrder(x1, x2)
            [startPos, finishPos] = calcOrder(x1,x2);
            % this calaculates the penalty of the line using the function normalLine with the parameters 
            % startPos, finishPos, opt, intersection, map, gradient, c
            intersection = normalLine(startPos, finishPos, opt, intersection, map, gradient, c);
        end
        % if the absolute value of x1 + x2 is greater than the absolute value of y1 and y2 and check is zero
        % and checkNaN is zero then this will use the go along the y axis to calaculate the x value
        if abs(x1 + x2) > abs(y1 + y2) & check == 0 & checkNaN == 0
            % sets opt to 1
            opt = 1;
            % calaculate the startPos and finishPos of the line using calcOrder(y1, y2)
            [startPos, finishPos] = calcOrder(y1,y2);
            % this calaculates the penalty of the line using the function normalLine with the parameters 
            % startPos, finishPos, opt, intersection, map, gradient, c
            intersection = normalLine(startPos, finishPos, opt, intersection, map, gradient, c);
        end
    end 
    % calaculates total error by adding pathDistance and intersection
    totalError = pathDistance + intersection;
    % sets fitness to 1 over totalError so the bigger the error the smaller the fitness
    fitness = 1/pathDistance + 1/intersection;
end

% This function checks everyCord on the line between two points and checks each point if they are on a 
% shape and adds a penalty if it is
function intersection = straightLine(startPos, finishPos, stationaryCoord, opt, intersection, map)
    % temp is equal to stationaryCoord
    temp = stationaryCoord;
    % this loops from startPos to finishPos
    for j = startPos:finishPos
        % movingPos is set to j
        movingPos = j;
        % opt equal 1 swap stationaryCoord and movingPos
        if opt == 1
            % stationaryCoord is set to movingPos
            stationaryCoord = movingPos;
            % movingPos is set to temp
            movingPos = temp;
        end
        % this checks if at position (stationaryCoord. movingPos) is equal to 0 (shape)
        if map(movingPos,stationaryCoord) == 0
            % this adds a penalty to line
            intersection = intersection + 2000;
        end
    end
end

% This function checks everyCord on the line between two points and checks each point if they are on a 
% shape and adds a penalty if it is
function intersection = normalLine(startPos, finishPos, opt, intersection, map, gradient, c)
    % this loops form startPos to finishPos
    for j = startPos:finishPos
        % sets xPos to j
        xPos = j;
        % calaculates the yPos using the equation of a line and then rounds the value
        yPos = round((gradient*j) + c);
        % if option is 1 then we calaculate the xPos
        if(opt == 1)
            % set yPos to j
            yPos = j;
            % calaculates the xPos doing the reverse of the equation of a line
            xPos = round((j - c) / gradient);
        end
        % this checks if at position (yPOs. xPos) is equal to 0 (shape)
        if map(yPos, xPos) == 0
            % this adds a penalty to line 
            intersection = intersection + 2000;
        end
    end
end

% This function requires two values and then it will return the smaller one as the startPos and the 
% bigger one as finishPos
function [startPos, finishPos] = calcOrder(value1, value2)
    % sets startPos to value1
    startPos = value1;
    % sets finishPos to value2
    finishPos = value2;
    if(value1 > value2)
        % if value1 is bigger than value2 then it will set startPos to value2
        startPos = value2;
        % sets finsihPos to value1
        finishPos = value1;
    end
end


% This function returns a index calaculated by doing a TournamentSelection with the number of rounds determined by k
function choice = TournamentSelection(k, population)
    % choice is set to -1
    choice = -1;
    % for 1 to k it will generate a random index from 1 to size of population 
    % Then it will check if choice == -1 or the fitness of the population at index ind
    % is greater than the fitness of the population at index choice. IF it is the set choice
    % equal to index
    for i=1:k
        ind = randi(size(population,1),1);
        if choice == -1 || population(ind, size(population,2)) > population(choice, size(population,i))
            choice = ind;
        end
    end
end

% Function for Bit Flip Mutation in a genetic algorithm
function temp_chromosome = BitFlipMutation(temp_chromosome)
    position1 = randi(size(temp_chromosome,2),1,1); % Generate a random position between 1 and 20
    position2 = randi(size(temp_chromosome,2),1,1); % Generate a random position between 1 and 20
    startPos = position1; % startPos is equal to position 1
    finsihPos = position2; % finishPos is equal to position 2
    if startPos > finsihPos
        startPos = position2; % Swap positions if startPos is greater than finsihPos
        finsihPos = position1;
    end
    distance = floor((finsihPos - startPos) / 2); % Calculate the distance between startPos and finsihPos
    count = 0; % Initialize a counter for the loop
    for i = startPos:finsihPos % Iterate over the range from startPos to finsihPos
        count = count + 1; % Increment the counter
        temp = temp_chromosome(i); % Store the value at position i
        temp_chromosome(i) = temp_chromosome(finsihPos - (i-startPos)); % Swap values using bit flip mutation
        temp_chromosome(finsihPos - (i-startPos)) = temp; % Swap values using bit flip mutation
        if distance == count
            break; % Exit the loop when the desired distance is reached
        end 
    end
end

function [temp_chromosome_1, temp_chromosome_2] = kpointCrossover(temp_chromosome_1, temp_chromosome_2, k)
    % Get the length of the chromosomes
    endPoint = size(temp_chromosome_1,2);

    % Initialize an array to store the crossover points
    crossOver = zeros(k,1);

    % Randomly select the first crossover point within the chromosome length
    crossOver(1) = randi(20, 1, 1);

    % Check if there are two crossover points
    if k == 2
        % Randomly select the second crossover point within the chromosome length
        crossOver(2) = randi(20, 1, 1);
        
        % Sort the crossover points in ascending order
        crossOver = sort(crossOver);
        
        % Update the endpoint for crossover
        endPoint = crossOver(2);
    end 

    % Initialize a temporary array
    temp = zeros(1,20);

    % Perform crossover between the selected points
    for i = crossOver:endPoint
        temp(i) = temp_chromosome_1(i);
        temp_chromosome_1(i) = temp_chromosome_2(i);
        temp_chromosome_2(i) = temp(i);
    end
end
function choice = RouletteWheelSelection(weights)
    % Calculate the cumulative sum of weights
    accumulation = cumsum(weights);
    % Generate a random number between 0 and 1
    p = rand();
    % Initialize the chosen index
    chosen_index = -1;
    % Iterate through the accumulation array to find the chosen index
    for index = 1 : length(accumulation)
        % if the cumulative sum at position index is bigger than p
        if (accumulation(index) > p)
            % Set the chosen index and exit the loop
            chosen_index = index;
            break;
        end
    end
    % Assign the chosen index as the output
    choice = chosen_index;
end

function [temp_chromosome_1, temp_chromosome_2] = UniformCrossover(temp_chromosome_1, temp_chromosome_2, swappingProb)
    % generate a random number from 0 to 1 20 times and store it in crossOver
    crossOver = rand(1,size(temp_chromosome_1,2));
    % creates a matrix of zeros which is 1 by 20
    temp = zeros(1,size(temp_chromosome_2,2));
    % iterates i from 1 to 20
    for i = 1:size(temp_chromosome_1,2)
        % if crossOver at position i is bigger than swappingProb
        if(crossOver(i) > swappingProb)
            % makes temp at position i equal to temp_chromosome_1 at position i
            temp(i) = temp_chromosome_1(i);
            % makes temp_chromosome_1 at position i equal to temp_chromosome_2 at position i
            temp_chromosome_1(i) = temp_chromosome_2(i);
            % makes temp_chromosome_2 at position i equal to temp at position i
            temp_chromosome_2(i) = temp(i);
        end
    end 
end

function temp_chromosome = ValueEncodingMutation(temp_chromosome,value)
    % generate a random number from 1 to 20 and store it in position1
    position = randi(size(temp_chromosome,2),1);
    % generates a random number between -value and value
    mutation_value = floor(rand() * 2 * value - value);
    while temp_chromosome(position) + mutation_value > 500 && temp_chromosome(position) + mutation_value <= 0
        mutation_value = floor(rand() * 2 * value - value);
        %mutation_value = floor(rand() * 2 * value - value);
    end
    temp_chromosome(position) = temp_chromosome(position) + mutation_value;
end

function selection = option1()
    error = 0;
    while error == 0
        try
            disp("Please choose a Selection");
            disp("0: Roulette wheel selection ");
            disp("1: Tournament selection ");
            disp("2: Rank-based Selection");
            disp("Please enter a number from 0 - 2 to choose the Selection");
            prompt = "Input:     ";
            selection = input(prompt, "s");
            selection = str2double(selection);
            error = 1;
        catch
            disp("Your Input was not a Valid Input");
            disp("Try again");
            error = 0;
        end
        if ~(isreal(selection) && rem(selection,1) == 0 && selection > -1 && selection < 3 && error ~= 0)
            disp("Your Input was not a Valid Input");
            disp("Try again");
            error = 0;
        end
    end
end 

function crossover = option2()
    error = 0;
    while error == 0
        try
            disp("Please choose a Crossover");
            disp("0: k-point Crossover ");
            disp("1: Uniform Crossover ");
            disp("Please enter a number from 0 - 1 to choose the Crossover");
            prompt = "Input:     ";
            crossover = input(prompt, "s");
            crossover = str2double(crossover);
            error = 1;
        catch
            disp("Your Input was not a Valid Input");
            disp("Try again");
            error = 0;
        end
        if ~(isreal(crossover) && rem(crossover,1) == 0 && crossover > -1 && crossover < 2 && error ~= 0)
            disp("Your Input was not a Valid Input");
            disp("Try again");
            error = 0;
        end
    end
end 

function mutation = option3()
    error = 0;
    while error == 0
        try
            disp("Please choose a Mutation");
            disp("0: bit-flip mutation");
            disp("1: Value Encoding");
            disp("Please enter a number from 0 - 1 to choose the Mutation");
            prompt = "Input:     ";
            mutation = input(prompt, "s");
            mutation = str2double(mutation);
            error = 1;
        catch
            disp("Your Input was not a Valid Input");
            disp("Try again");
            error = 0;
        end
        if ~(isreal(mutation) && rem(mutation,1) == 0 && mutation > -1 && mutation < 2 && error ~= 0)
            disp("Your Input was not a Valid Input");
            disp("Try again");
            error = 0;
        end
    end
end 
