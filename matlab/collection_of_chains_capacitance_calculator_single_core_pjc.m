numwalkers = 10000; % Number of random walkers to use

number_of_chains = number_of_chains_to_generate;

capacitances = zeros(1,number_of_chains); % Capacitances of chains

did_it_hit = zeros(number_of_chains,numwalkers); % Did each walker hit
chain_time = zeros(number_of_chains); % How long to do each chain
hit_record = zeros(numwalkers, 3, number_of_chains);
absorbing_layer = 0.0008;

L=0;

for chain = 1:number_of_chains
    
    hits = 0;
    losts = 0;
    
    % Use the smallest R that includes all points
    x_max = max(max(collection_of_chains(:,1,:,chain)));
    x_min = min(min(collection_of_chains(:,1,:,chain)));
    y_max = max(max(collection_of_chains(:,2,:,chain)));
    y_min = min(min(collection_of_chains(:,2,:,chain)));
    z_max = max(max(collection_of_chains(:,3,:,chain)));
    z_min = min(min(collection_of_chains(:,3,:,chain)));
    x_mags = [sqrt(x_max*x_max), sqrt(x_min*x_min)];
    y_mags = [sqrt(y_max*y_max), sqrt(y_min*y_min)];
    z_mags = [sqrt(z_max*z_max), sqrt(z_min*z_min)];
    x_max_mag = max(x_mags);
    y_max_mag = max(y_mags);
    z_max_mag = max(z_mags);
    R = sqrt(x_max_mag*x_max_mag + y_max_mag*y_max_mag + z_max_mag*z_max_mag);
    rinfinity = 500*R;
        
    tic
    
    % Use parallel for loop to run walkers simultaneously
    for walkers = 1:numwalkers
        shortstepsize = 0.01;
        distances = zeros(1, number_of_plates_in_chain);
        
        lost = false;
        hit = false;
        
        % Choose random position on sphere R to begin
        u1=rand;
        u2=rand;
        theta=2*3.142*u1;
        phi=acos(1-2*u2);
        x=R*cos(theta)*sin(phi);
        y=R*sin(theta)*sin(phi);
        z=R*cos(phi);
        r=R;
        
        % Loop until walker either hits crystal or escapes past rinfinity
        while (lost==false)&(hit==false)
        
            if r > R
                dr = r - R;
            else
                dr = shortstepsize;
            end
            
            % Move by distance dr in random direction
            u1=rand;
            u2=rand;
            theta=2*3.142*u1;
            phi=acos(1-2*u2);
            dx=dr*cos(theta)*sin(phi);
            dy=dr*sin(theta)*sin(phi);
            dz=dr*cos(phi);
            x=x+dx; y=y+dy; z=z+dz;
            r=sqrt(x*x+y*y+z*z);
        
        
            if r > rinfinity
                lost = true;
                losts=losts+1;
                continue
            elseif r > R
                continue
            else
                
                % Test whether walker is within a distance absorbing_layer of surface of
                % surface of any of the plates in the chain
                %If it is, hit = true, exit loop
                for current_plate = 1:number_of_plates_in_chain
                      p = [x y z];
                      shape_points = zeros(12,3);
                      q = p - centres(chain,:,current_plate);
                      for i = 1:12
                          for j = 1:3
                        shape_points(i,j) = collection_of_chains(i,j,current_plate,chain) - centres(chain,j,current_plate);
                          end
                      end
                      %[in, distances(current_plate)] = new_gjk(p, collection_of_chains(:,:,current_plate,chain));
                      [in, distances(current_plate)] = new_gjk(q, shape_points);

                      if in == true || distances(current_plate) < absorbing_layer
                          if in==true
                              in = in
                          end
                          hit = true;
                          did_it_hit(chain,walkers) = 1;
                          hits = hits+1;
                          hit_record(walkers,1,chain) = x;
                          hit_record(walkers,2,chain) = y;
                          hit_record(walkers,3,chain) = z;
                          break
                      end
                       
                end
            end
            
            % If it didn't hit and is inside R, step size is the distance
            % to the surface
            shortstepsize = min(distances);
            
        end
        
    chain
     walkers
     ratio = hits/(hits+losts)
     capacitance = (hits / (hits + losts))*R
    end
    
    toc
    chain_time(chain) = toc;
    
     
    % capac_record(walkers,chain) = (hits / (hits + losts))*R;
    
    % After all walkers complete, add up capacitance for chain
   % for i=1:numwalkers
     %   capacitances(chain) = capacitances(chain) + R*did_it_hit(chain,i)/numwalkers;
   % end
end

