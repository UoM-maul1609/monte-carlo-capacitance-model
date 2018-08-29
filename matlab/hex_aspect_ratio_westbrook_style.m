numwalkers = 10000;

number_of_chains = number_of_chains_to_generate;

%hit_record = zeros(numwalkers, 3, number_of_chains);
%capac_record = zeros(numwalkers, 1, number_of_chains);
capacitances = zeros(1,number_of_chains);

did_it_hit = zeros(number_of_chains,numwalkers);
chain_time = zeros(number_of_chains);

L=0;

for chain = 1:number_of_chains
    
    L= L+chain*0.1;
    
    R = sqrt(a*a + L*L);
    %R = 1.5;
    rinfinity = 500*R;
    
    hits = 0;
    losts = 0;

    numsofar = 0;
    
    x_max = max(max(collection_of_chains(:,1,:,chain)));
    x_min = min(min(collection_of_chains(:,1,:,chain)));
    
    y_max = max(max(collection_of_chains(:,2,:,chain)));
    y_min = min(min(collection_of_chains(:,2,:,chain)));
    
    z_max = max(max(collection_of_chains(:,3,:,chain)));
    z_min = min(min(collection_of_chains(:,3,:,chain)));

    tic
    
    parfor walkers = 1:numwalkers
        shortstepsize = 0.01;
        distances = zeros(1, number_of_plates_in_chain);
        
        lost = false;
        hit = false;
    
        u1=rand;
        u2=rand;
        theta=2*3.142*u1;
        phi=acos(1-2*u2);
        x=R*cos(theta)*sin(phi);
        y=R*sin(theta)*sin(phi);
        z=R*cos(phi);
        r=R;
    
        while (lost==false)&(hit==false)
        
            if r > R
                dr = r - R;
            else
                dr = shortstepsize;
            end
        
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
                %losts = losts + 1;
                continue
            elseif r > R
                continue
            else
                for current_plate = 1:number_of_plates_in_chain
                      p = [x y z];
                      
                      %[in, distances(current_plate)] = new_gjk(p, collection_of_chains(:,:,current_plate,chain));
                      distances(current_plate) = real_distance_1_final(p, collection_of_chains(:,:,current_plate,chain));
                      %shortstepsize = 0.01;
                     % is_in = in
                      current_distance = distances(current_plate)
                      if min(distances) <= 0.001
                          hit = true;
                          did_it_hit(chain,walkers) = 1;
                          %hits = hits+1;
                          %hit_record(walkers,1,chain) = x;
                          %hit_record(walkers,2,chain) = y;
                          %hit_record(walkers,3,chain) = z;
                          break 
                      end  
                end
            end
            
            shortstepsize = min(distances)+0.01
            
        end
        
        

        %chain
        %walkers
        %ratio = hits/(hits+losts)
        %capacitance = (hits / (hits + losts))*R
        %capac_record(walkers,chain) = (hits / (hits + losts))*R;
    
    end
    
    toc
    chain_time(chain) = toc;
    
    for i=1:numwalkers
        capacitances(chain) = capacitances(chain) + R*did_it_hit(chain,i)/numwalkers;
    end
    %capacitances(chain) = (hits / (hits + losts))*R;
end
