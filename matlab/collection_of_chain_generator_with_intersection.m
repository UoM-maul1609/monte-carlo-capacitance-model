a = 1;
L = 0.5;
number_of_plates_in_chain = 5;

number_of_chains_to_generate = 1;

collection_of_chains = zeros(12, 3, number_of_plates_in_chain, number_of_chains_to_generate);
aggregation_indices = zeros(1,number_of_chains_to_generate);

x_axis_extension_bias = 0;

%collection_of_chains = zeros(12, 3, number_of_plates_in_chain);
centres = zeros(number_of_chains_to_generate,3,number_of_plates_in_chain);
D_maxes = zeros(1,number_of_chains_to_generate);
basic_hex = [-a/2 -a*sqrt(3)/2 L/2; a/2 -a*sqrt(3)/2 L/2; a 0 L/2; a/2 a*sqrt(3)/2 L/2; -a/2 a*sqrt(3)/2 L/2;  -a 0 L/2; -a/2 -a*sqrt(3)/2 -L/2; a/2 -a*sqrt(3)/2 -L/2; a 0 -L/2; a/2 a*sqrt(3)/2 -L/2; -a/2 a*sqrt(3)/2 -L/2;  -a 0 -L/2];

for current_chain = 1:number_of_chains_to_generate
    tic
collection_of_chains(:,:,1, current_chain) = basic_hex;

for hexnum = 2:number_of_plates_in_chain
    
    %Choose random orientation of new hex plate
    theta_x = 2*3.142*rand;
    theta_y = 2*3.142*rand;
    theta_z = 2*3.142*rand;
    collection_of_chains(:,:,hexnum,current_chain) = rotate_shape_rad(basic_hex, theta_x, theta_y, theta_z);
     % Move new hex to centre of previous hexagon
    collection_of_chains(:,:,hexnum,current_chain) = translate_shape(collection_of_chains(:,:,hexnum,current_chain), centres(1,1,hexnum-1), centres(1,2,hexnum-1), centres(1,3,hexnum-1));
    centres(current_chain,:,hexnum) = translate_shape(centres(current_chain,:,hexnum), x_axis_extension_bias, 0 , 0);
    % Move new plate in x direction by bias amount
    collection_of_chains(:,:,hexnum,current_chain) = translate_shape(collection_of_chains(:,:,hexnum,current_chain), x_axis_extension_bias, 0, 0);
   
    % Choose random theta, phi for movement of hexagon
    translation_theta = 2*3.142*rand;
    translation_phi = acos(1-2*rand);
    translation_length = 0.05;
    % Array records if it intersects any of the previous hexagons in chain
    in_chain = zeros(1,hexnum-1);
    % Array records if it intersected in this or the last iteration
    in_iteration = zeros(1,2);
    
    % First iteration: 
    % Translate one step in chosen direction
    translation_x = translation_length*sin(translation_theta)*cos(translation_phi);
    translation_y = translation_length*sin(translation_theta)*sin(translation_phi);
    translation_z = translation_length*cos(translation_theta);
    collection_of_chains(:,:,hexnum,current_chain) = translate_shape(collection_of_chains(:,:,hexnum,current_chain), translation_x, translation_y, translation_z);
    centres(current_chain,:,hexnum) = translate_shape(centres(current_chain,:,hexnum-1), translation_x, translation_y, translation_z);
    % in_iteration(2) records how many of the pre-existing hexagons it
    % intersects, if zero, reverse direction
    for hex_to_test = 1:hexnum-1
        in_chain(hex_to_test) = new_gjk_polyhedra(collection_of_chains(:,:,hexnum,current_chain), collection_of_chains(:,:,hex_to_test,current_chain));
    end
    for m = 1:hexnum-1
        in_iteration(2) = in_iteration(2) + in_chain(m);
    end
    
    if(in_iteration(2) == 0)
        translation_length = -translation_length;
    end
    correct = false;
    while(correct == false)
        
        in_iteration(1) = in_iteration(2);
        in_iteration(2) = 0;
        
        translation_x = translation_length*sin(translation_theta)*cos(translation_phi);
        translation_y = translation_length*sin(translation_theta)*sin(translation_phi);
        translation_z = translation_length*cos(translation_theta);
        collection_of_chains(:,:,hexnum,current_chain) = translate_shape(collection_of_chains(:,:,hexnum,current_chain), translation_x, translation_y, translation_z);
        centres(current_chain,:,hexnum) = translate_shape(centres(current_chain,:,hexnum-1), translation_x, translation_y, translation_z);
        
        for hex_to_test = 1:hexnum-1
            in_chain(hex_to_test) = new_gjk_polyhedra(collection_of_chains(:,:,hexnum,current_chain), collection_of_chains(:,:,hex_to_test,current_chain));
        end
        
        for n = 1:hexnum-1
            in_iteration(2) = in_iteration(2) + in_chain(n);
        end
        
        if( (in_iteration(2) == 0) && (in_iteration(1) > 0) )
            translation_length = -2*translation_length;
            translation_x = translation_length*sin(translation_theta)*cos(translation_phi);
            translation_y = translation_length*sin(translation_theta)*sin(translation_phi);
            translation_z = translation_length*cos(translation_theta);
            collection_of_chains(:,:,hexnum,current_chain) = translate_shape(collection_of_chains(:,:,hexnum,current_chain), translation_x, translation_y, translation_z);
            centres(current_chain,:,hexnum) = translate_shape(centres(current_chain,:,hexnum-1), translation_x, translation_y, translation_z);
            correct = true;
        end
        
    end
end

% Finds current centre of mass of aggregate, then translates all hexagons
% and points back to be centred around the origin
centre_of_mass = (1/number_of_plates_in_chain)*sum(centres, 3);

collection_of_chains(:,1,:,current_chain) = collection_of_chains(:,1,:,current_chain) - centre_of_mass(1);
collection_of_chains(:,2,:,current_chain) = collection_of_chains(:,2,:,current_chain) - centre_of_mass(2);
collection_of_chains(:,3,:,current_chain) = collection_of_chains(:,3,:,current_chain) - centre_of_mass(3);

all_distances = zeros(12,number_of_plates_in_chain,12,number_of_plates_in_chain);
for i = 1:12
    for j = 1:number_of_plates_in_chain
        for k = 1:12
            for l = 1:number_of_plates_in_chain
                all_distances(i,j,k,l) = sqrt( (collection_of_chains(i,1,j,current_chain) - collection_of_chains(k,1,l,current_chain))^2 + (collection_of_chains(i,2,j,current_chain) - collection_of_chains(k,2,l,current_chain))^2 + (collection_of_chains(i,3,j,current_chain) - collection_of_chains(k,3,l,current_chain))^2);
            end
        end
    end
end
D_maxes(current_chain) = max(max(max(max(all_distances))));

rect_face_definitions = [1 2 8 7; 2 3 9 8; 3 4 10 9; 4 5 11 10; 5 6 12 11; 6 1 7 12];
hex_face_definitions = [1 2 3 4 5 6; 7 8 9 10 11 12];

distance_sum_for_aggregate = 0;
max_distance_sum_for_aggregate = 0;
for i = 1:number_of_plates_in_chain
    for j = 1:number_of_plates_in_chain
        max_distance_sum_for_aggregate = max_distance_sum_for_aggregate + abs(i-j)*2*a;
    end
end
for i = 1:number_of_plates_in_chain
    for j = 1:number_of_plates_in_chain
       distance_sum_for_aggregate = distance_sum_for_aggregate + sqrt( (centres(current_chain,1,i) - centres(current_chain,1,j))^2 + (centres(current_chain,2,i) - centres(current_chain,2,j))^2 + (centres(current_chain,3,i) - centres(current_chain,3,j))^2);
    end
end

aggregation_indices(current_chain) = distance_sum_for_aggregate/max_distance_sum_for_aggregate;

toc

end