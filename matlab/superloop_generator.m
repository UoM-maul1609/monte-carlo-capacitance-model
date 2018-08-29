a = 1;
L = 0.5;

%number_of_plates_in_chain = ;

number_of_chains_to_generate = 65;
D_maxes = zeros(1,number_of_chains_to_generate);
% Array storing all chains generated
collection_of_chains = zeros(12, 3, number_of_plates_in_chain, number_of_chains_to_generate);
aggregation_indices = zeros(1,number_of_chains_to_generate);
%Array storing positions of centres of hexagons
centres = zeros(number_of_chains_to_generate,3,number_of_plates_in_chain);

% Basic hexagonal plate centred on origin from which others are generated
basic_hex = [L/2 -a*sqrt(3)/2 -a/2; L/2 -a*sqrt(3)/2 a/2; L/2 0 a; L/2 a*sqrt(3)/2 a/2; L/2 a*sqrt(3)/2 -a/2; L/2 0 -a; -L/2 -a*sqrt(3)/2 -a/2; -L/2 -a*sqrt(3)/2 a/2; -L/2 0 a; -L/2 a*sqrt(3)/2 a/2; -L/2 a*sqrt(3)/2 -a/2; -L/2 0 -a; ];

%Parameters that set how isotropic the distributions of randomly chosen
%parameters are. (alpha = beta in beta distribution). 1 = uniformly random,
% -> 0 = isotropic distribution.

beta_x = 1; % Controls orientation rotation about x-axis
beta_y = 1; % Controls orientation rotation about y-axis
beta_z = 1; % Controls orientation rotation about z-axis
beta_trans = 1; % Controls translation angle around z-axis


% Only varies beta_trans until reaches 0.7
for current_chain = 1:15
    
    tic % Starts timer
    
    collection_of_chains(:,:,1, current_chain) = basic_hex; % First in chain is basic hex

    for current_hex = 2:number_of_plates_in_chain

        % Choose random orientation of new hex plate
        % Rotate first around x-axis, then y, then z
        % Angles chosen from beta-distribution with parameters given above.
        % Beta distribution gives control of isotropy
        theta_x = 2*3.142*betarnd(beta_x,beta_x);
        theta_y = 2*3.142*betarnd(beta_y,beta_y);
        theta_z = 2*3.142*betarnd(beta_z,beta_z);
        collection_of_chains(:,:,current_hex,current_chain) = rotate_shape_rad(basic_hex, theta_x, theta_y, theta_z);
        
        % Choose random theta, phi for movement of hexagon
        % phi uniformly distributed, theta chosen from beta distribution
        translation_phi = 2*3.142*rand;
        translation_theta = 3.142*betarnd(beta_trans,beta_trans);
        
             
        % Moves hex by translation_length in each step, then tests for
        % intersection
        translation_length = 0.05;
        % Array records if it intersects any of the previous hexagons in chain
        in_chain = zeros(1,current_hex-1);
        % Array records if it intersected in this or the last iteration
        in_iteration = zeros(1,2);

        % First iteration: 
        % Translate one step in chosen direction
        translation_x = translation_length*sin(translation_theta)*cos(translation_phi);
        translation_y = translation_length*sin(translation_theta)*sin(translation_phi);
        translation_z = translation_length*cos(translation_theta);
        collection_of_chains(:,:,current_hex,current_chain) = translate_shape(collection_of_chains(:,:,current_hex,current_chain), translation_x, translation_y, translation_z);
        centres(current_chain,:,current_hex) = translate_shape(centres(current_chain,:,current_hex), translation_x, translation_y, translation_z);
        % in_iteration(2) records how many of the pre-existing hexagons it
        % intersects, if zero, reverse direction
        for hex_to_test = 1:current_hex-1
            in_chain(hex_to_test) = new_gjk_polyhedra(collection_of_chains(:,:,current_hex,current_chain), collection_of_chains(:,:,hex_to_test,current_chain));
        end
        for m = 1:current_hex-1
            in_iteration(2) = in_iteration(2) + in_chain(m);
        end
        
        % If doesn't intersect after first movement, move in other
        % direction
        if(in_iteration(2) == 0)
            translation_length = -translation_length;
        end
        
        % Loop until hex's just intersect
        correct = false;
        while(correct == false)

            in_iteration(1) = in_iteration(2);
            in_iteration(2) = 0;

            translation_x = translation_length*sin(translation_theta)*cos(translation_phi);
            translation_y = translation_length*sin(translation_theta)*sin(translation_phi);
            translation_z = translation_length*cos(translation_theta);
            collection_of_chains(:,:,current_hex,current_chain) = translate_shape(collection_of_chains(:,:,current_hex,current_chain), translation_x, translation_y, translation_z);
            centres(current_chain,:,current_hex) = translate_shape(centres(current_chain,:,current_hex), translation_x, translation_y, translation_z);
            
            % Check if hex intersects any of the previous hex's
            for hex_to_test = 1:current_hex-1
                in_chain(hex_to_test) = new_gjk_polyhedra(collection_of_chains(:,:,current_hex,current_chain), collection_of_chains(:,:,hex_to_test,current_chain));
            end
            
            % Add up how many it intersects
            for n = 1:current_hex-1
                in_iteration(2) = in_iteration(2) + in_chain(n);
            end
            
            % If hex's don't intersect and did on last iteration, move hex
            % back two steps so they overlap a bit, then we are done
            if( (in_iteration(2) == 0) && (in_iteration(1) > 0) )
                translation_length = -2*translation_length;
                translation_x = translation_length*sin(translation_theta)*cos(translation_phi);
                translation_y = translation_length*sin(translation_theta)*sin(translation_phi);
                translation_z = translation_length*cos(translation_theta);
                collection_of_chains(:,:,current_hex,current_chain) = translate_shape(collection_of_chains(:,:,current_hex,current_chain), translation_x, translation_y, translation_z);
                centres(current_chain,:,current_hex) = translate_shape(centres(current_chain,:,current_hex), translation_x, translation_y, translation_z);
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

    centres(current_chain,1,:) = centres(current_chain,1,:) - centre_of_mass(1);
    centres(current_chain,2,:) = centres(current_chain,2,:) - centre_of_mass(2);
    centres(current_chain,3,:) = centres(current_chain,3,:) - centre_of_mass(3);
      
    
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
    
    % aggregation index = SUM_OVER_I_AND_J[ (distance between centres i and j) / (max
    % possible distance between centres i and j) ]
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

    beta_trans = beta_trans - 0.01;
end

for current_chain = 16:number_of_chains_to_generate
    
    tic % Starts timer
    
    collection_of_chains(:,:,1, current_chain) = basic_hex; % First in chain is basic hex

    for current_hex = 2:number_of_plates_in_chain

        % Choose random orientation of new hex plate
        % Rotate first around x-axis, then y, then z
        % Angles chosen from beta-distribution with parameters given above.
        % Beta distribution gives control of isotropy
        theta_x = 2*3.142*betarnd(beta_x,beta_x);
        theta_y = 2*3.142*betarnd(beta_y,beta_y);
        theta_z = 2*3.142*betarnd(beta_z,beta_z);
        collection_of_chains(:,:,current_hex,current_chain) = rotate_shape_rad(basic_hex, theta_x, theta_y, theta_z);
        
        % Choose random theta, phi for movement of hexagon
        % phi uniformly distributed, theta chosen from beta distribution
        translation_phi = 2*3.142*rand;
        translation_theta = 3.142*betarnd(beta_trans,beta_trans);
        
             
        % Moves hex by translation_length in each step, then tests for
        % intersection
        translation_length = 0.05;
        % Array records if it intersects any of the previous hexagons in chain
        in_chain = zeros(1,current_hex-1);
        % Array records if it intersected in this or the last iteration
        in_iteration = zeros(1,2);

        % First iteration: 
        % Translate one step in chosen direction
        translation_x = translation_length*sin(translation_theta)*cos(translation_phi);
        translation_y = translation_length*sin(translation_theta)*sin(translation_phi);
        translation_z = translation_length*cos(translation_theta);
        collection_of_chains(:,:,current_hex,current_chain) = translate_shape(collection_of_chains(:,:,current_hex,current_chain), translation_x, translation_y, translation_z);
        centres(current_chain,:,current_hex) = translate_shape(centres(current_chain,:,current_hex), translation_x, translation_y, translation_z);
        % in_iteration(2) records how many of the pre-existing hexagons it
        % intersects, if zero, reverse direction
        for hex_to_test = 1:current_hex-1
            in_chain(hex_to_test) = new_gjk_polyhedra(collection_of_chains(:,:,current_hex,current_chain), collection_of_chains(:,:,hex_to_test,current_chain));
        end
        for m = 1:current_hex-1
            in_iteration(2) = in_iteration(2) + in_chain(m);
        end
        
        % If doesn't intersect after first movement, move in other
        % direction
        if(in_iteration(2) == 0)
            translation_length = -translation_length;
        end
        
        % Loop until hex's just intersect
        correct = false;
        while(correct == false)

            in_iteration(1) = in_iteration(2);
            in_iteration(2) = 0;

            translation_x = translation_length*sin(translation_theta)*cos(translation_phi);
            translation_y = translation_length*sin(translation_theta)*sin(translation_phi);
            translation_z = translation_length*cos(translation_theta);
            collection_of_chains(:,:,current_hex,current_chain) = translate_shape(collection_of_chains(:,:,current_hex,current_chain), translation_x, translation_y, translation_z);
            centres(current_chain,:,current_hex) = translate_shape(centres(current_chain,:,current_hex), translation_x, translation_y, translation_z);
            
            % Check if hex intersects any of the previous hex's
            for hex_to_test = 1:current_hex-1
                in_chain(hex_to_test) = new_gjk_polyhedra(collection_of_chains(:,:,current_hex,current_chain), collection_of_chains(:,:,hex_to_test,current_chain));
            end
            
            % Add up how many it intersects
            for n = 1:current_hex-1
                in_iteration(2) = in_iteration(2) + in_chain(n);
            end
            
            % If hex's don't intersect and did on last iteration, move hex
            % back two steps so they overlap a bit, then we are done
            if( (in_iteration(2) == 0) && (in_iteration(1) > 0) )
                translation_length = -2*translation_length;
                translation_x = translation_length*sin(translation_theta)*cos(translation_phi);
                translation_y = translation_length*sin(translation_theta)*sin(translation_phi);
                translation_z = translation_length*cos(translation_theta);
                collection_of_chains(:,:,current_hex,current_chain) = translate_shape(collection_of_chains(:,:,current_hex,current_chain), translation_x, translation_y, translation_z);
                centres(current_chain,:,current_hex) = translate_shape(centres(current_chain,:,current_hex), translation_x, translation_y, translation_z);
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

    
    centres(current_chain,1,:) = centres(current_chain,1,:) - centre_of_mass(1);
    centres(current_chain,2,:) = centres(current_chain,2,:) - centre_of_mass(2);
    centres(current_chain,3,:) = centres(current_chain,3,:) - centre_of_mass(3);
    
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
    
    % aggregation index = SUM_OVER_I_AND_J[ (distance between centres i and j) / (max
    % possible distance between centres i and j) ]
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

    beta_trans = beta_trans - 0.002;
    beta_x = beta_x - 0.005;
    beta_y = beta_y - 0.005;
end