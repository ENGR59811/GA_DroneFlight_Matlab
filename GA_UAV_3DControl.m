%% Written by Kelvin Ma                                                                     
% City College of New York, CUNY                                                                            
% Date: August 16, 2020
% UAV Control with Genetic Algorithm 
% Given the number of UAVs in a swarm and a radius of communication, 
% use a genetic algorithm to uniformly
% spread the swarm such that they 
% retain connectivity with its neighbors while maintaining a distance  
% close to radius from all its neighbors. 
clc; clear;
% User inputs and initalizations
% 
%user inputs:
noof_UAVs = 10; % number of UAVs = 20
% good: speed = 0.1* rcom: rcom = 10, speed = 1
% bad: when speed == rcom: rcom = 10, speed = 10
% ugly: speed = 0.8* rcom: rcom = 10, speed = 8
speed = 1; % constant speed of each UAV
rcom = 10; % radius of communication 
noof_moves = 10; % 20 total moves for each UAV to make 
% Create a 3D matrix to store all UAV positions - initalized all to 0's
% UAV_positions[noof_UAVs][3][noof_moves] 
% 3 is for spatial (x, y, z) coordinates
UAV_positions = zeros(noof_UAVs, 3, noof_moves); 

UAV_positions(:,:,1) = (rcom/sqrt(2)).*rand(noof_UAVs, 3);

% Give each UAV a random starting location within half rcom of 0
% rand function creates a random matrix of dimensions rand_mat[noof_UAVs][3]
% where each element is between 0 and 1, so we multiply by rcom/sqrt(2) 
% so that the distance away from (0,0,rcom) is always less than rcom.
% UAV_positions is a 3-d matrix: [uav_no][x,y,z][time]
% Example: for positions of UAVs at time 1:
%     5.7610    1.1145   14.6368 <- UAV_1
%     6.4049    6.8631   10.2525 <- UAV_2
%     0.8979    6.7682   16.0043
%     6.4585    3.4321   16.6043
%     4.4715    5.6588   14.7994
%     0.6897    1.0033   15.3580
%     1.9693    2.9823   15.2547
%     3.8670    6.4752   12.7735
%     6.7706    5.6018   14.6349
%     6.8228    6.7846   11.2105 <- UAV_10

% we will assume that we want the UAVs to begin GA at a altitude higher
% than rcom (This makes it less likely to crash in the ground)
UAV_positions(:,3,1) = UAV_positions(:,3,1)+rcom;
% Chromosomes
% The chromosome will represent the direction for the UAV to move
% Each chromosome will be represented as a vector of length 5 
% Where the most significant 2 bits will represent altitude as such:
% [0, 0, x, x, x] = Stationary 
% [0, 1, x, x, x] = Down 
% [1, 0, x, x, x] = Up
% [1, 1, x, x, x] = Maintain current altitude
% There are 4 cardinal directions and 4 ordinal directions - they are:
% North, Northeast, East, Southeast, South, Southwest, West, Northwest
% we can encode with the least significant 3 bits of the chromosome:
%   [x, x, 0,0,0] = North   [x, x, 0,0,1] = North-East
%   [x, x, 0,1,0] = East    [x, x, 0,1,1] = South-East
%   [x, x, 1,0,0] = South   [x, x, 1,0,1] = South-West
%   [x, x, 1,1,0] = West    [x, x, 1,1,1] = North-West
%
% Since there are no restrictions on the chromosomes other than length,
% there is no need for the A and b constraint matrices.
% The fitness function for these chromosomes are at the end of the script
chromosome_length = 5;
Lb = zeros(1,chromosome_length); % lower bound of gene values
Ub = ones(1,chromosome_length); % upper bound of gene values 
% indices of the chromosome which are integers (all of them in this case)
int_indices = 1:chromosome_length; 
% Displacement vectors on the X-Y plane 
% The chromosomes are the DIRECTION for the UAV to move. Each direction and
% speed generate a displacement vector in space for the UAV. The
% magnitude of this displacement vector is always equal to speed. 
%
% Using the encodings for each direction, we can create an 2D array 
% that holds the unit displacement vectors, then multiply all the unit
% vectors by speed. The array dimensions are 
% displacement_vector[8][2] - 8 possible directions, 2 for x displacement
% and y displacement.
%
% for example, using the three least significant bits:
% (index) [unit displacement vector] 
% 1=(0+1)  [0, 1]                       - north      (xx000 is binary 0)
% 2=(1+1)  [1/sqrt(2), 1/sqrt(2)]       - northeast  (xx001 is binary 1)
% 3=(2+1)  [1, 0]                       - east       (xx010 is binary 2)
% 4=(3+1)  [1/sqrt(2), -1/sqrt(2)]      - southeast  (xx011 is binary 3)
% 5=(4+1)  [0, -1]                      - south      (xx100 is binary 4)
% 6=(5+1)  [-1/sqrt(2)), -1/sqrt(2)]    - southwest  (xx101 is binary 5) 
% 7=(6+1)  [-1, 0]                      - west       (xx110 is binary 6)
% 8=(7+1)  [-1/sqrt(2)), 1/sqrt(2)]     - northwest  (xx110 is binary 7) 
% Note: the +1 is because MATLAB vectors start with index 1
% chromosome : b2b1b0 -> convert to decimal 0-7 and use as an index (+1):
% [N, NE, E, SE, S, SW, W, NW]
% chromosome 101 -> 5 -> 5+1 -> SW
displacement_vectors = [[0, 1]; [1/sqrt(2), 1/sqrt(2)]; ...
                        [1, 0]; [1/sqrt(2), -1/sqrt(2)]; ... 
                        [0, -1]; [-1/sqrt(2), -1/sqrt(2)]; ... 
                        [-1, 0]; [-1/sqrt(2), 1/sqrt(2)]]; 
displacement_vectors = speed .* displacement_vectors;
% Running the Genetic Algorithm for all UAV from time 0 to end
% set options to: do not display (change off to iter to display steps) 
options = optimoptions('ga','display','off');  % set not to display

% Each UAV must run their own genetic algorithms at every time step
for t = 1:noof_moves-1
    fprintf("\n\n***** TIME IS NOW: %d*****\n", t)
    % iterate through the UAVs, the fitness function will change depending
    % on the UAV number we are on
    for UAV_num = 1:noof_UAVs
        fprintf("***** BEGINNING GA FOR UAV #%d *****\n", UAV_num)
        % get the neighbors of the UAV that is running GA - exclude itself 
        neighbors = UAV_positions(:,:,t);
        neighbors(UAV_num,:) = [];
        % get current position of UAV that is running ga
        current_position = UAV_positions(UAV_num, :, t); 
        
        % create a function pointer for the fitness function
        % (in Matlab we declare this inside the loop, o.w. does not update
        % parameter values at each loop iteration)
        fit_func = @(chromosome) fitness_function(chromosome, ...
                   current_position, speed, neighbors, rcom, ...
                   displacement_vectors);
        
        % Run the genetic algorithm for the current node
        % Matlab has its own default values for population size and number
        % of generations
        % selection contains the best chromosome returned by ga (it is the
        % next direction to move):
        selection = ga(fit_func,chromosome_length,[],[], ...
                    [],[],Lb,Ub,[],int_indices, options);
                
        % calculate the next position based on selection from ga:
        
        displacement=displacement_vectors(bin2dec(strrep(num2str(selection(3:end)), ' ', ''))+1,:);
        next_position = get_next_position(selection, current_position, ... 
                                            speed, displacement_vectors);
        
        % keep a record of these positions for the next ga
        UAV_positions(UAV_num,:,t+1) = next_position;
        fprintf("\tselection fitness = %f\n", fit_func(selection))
        fprintf("***** FINISHED GA FOR UAV #%d AT TIME %d *****\n", ...
            UAV_num, t)
    end
end

% Visualize the UAVs
visualize_UAVs3D(UAV_positions, noof_moves, noof_UAVs, rcom);
open_anim = true;
while open_anim == true
    fprintf("Replay the animation?\n\tenter 1 for yes\n\tenter 2 for no\n");
    user_choice = input('');
    if user_choice == 2
        open_anim = false;
        fprintf('Thank you. Have a nice day!\n');
        break;
    else
        close all
        visualize_UAVs3D(UAV_positions, noof_moves, noof_UAVs, rcom);
    end
end
% Fitness function
% Functions in MATLAB which are multiple lines are defined at the bottom
% of the script. 

% MATLAB will use GA to minimize the value variable called fitness_score
function fitness_score = fitness_function(chromosome, position, speed, ...
                            neighbors, rcom, displacement_vectors)
    fitness_score = 0;
    neighbor_count = 0;
    % calculate the candidate position from the chromosome and displacement 
    % vectors
    candidate_pos = get_next_position(chromosome, position, speed, ...
        displacement_vectors);
    
    % iterate through all the neighbors and calculate the fitness
    % fitness += (rcom - distance from neighbor[i])
    % ga minimizes this number, the closer to 0 the better 
    for i = 1:length(neighbors)
       distance = norm(candidate_pos - neighbors(i,:)); 
       if distance <= rcom
           neighbor_count = neighbor_count + 1;
           fitness_score = fitness_score + (rcom - distance);
       end
    end
    % maximum penalty for positions with no neighbors or positions that
    % crash to the ground
    if neighbor_count < 1 || candidate_pos(3) <= rcom
        fitness_score = abs(intmax) ;
    end 
end
% Getting the next position from the chromosome
% Functions in MATLAB which are multiple lines are defined at the bottom
% of the script. 
% this function gets the next position from a chromosome 
function next_pos = get_next_position(chromosome, position, speed, ...
                                        displacement_vectors)


   
    if bin2dec(strrep(num2str(chromosome(1:2)), ' ', '')) == 0 
        % stationary - next position = current_position 
        next_pos = position;
    elseif bin2dec(strrep(num2str(chromosome(1:2)), ' ', '')) == 1 
        % UAV moves up 
        
        displacement = displacement_vectors(bin2dec(strrep(num2str(chromosome(3:end)), ' ', ''))+1,:);
        next_pos(1:2)= position(1:2) + displacement;
        next_pos(3) = position(3) + speed; 
    elseif bin2dec(strrep(num2str(chromosome(1:2)), ' ', '')) == 2
        % UAV moves down 
        displacement = displacement_vectors(bin2dec(strrep(num2str(chromosome(3:end)), ' ', ''))+1,:);
        next_pos(1:2)= position(1:2) + displacement;
        next_pos(3) = position(3) - speed; 
    elseif bin2dec(strrep(num2str(chromosome(1:2)), ' ', '')) == 3 
        % UAV mainstains altitude 
        next_pos = position;
        displacement = displacement_vectors(bin2dec(strrep(num2str(chromosome(3:end)), ' ', ''))+1,:);
        next_pos(1:2)= position(1:2) + displacement;
    else
    end
end