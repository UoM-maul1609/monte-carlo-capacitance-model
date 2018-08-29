numwalkers = 100000;

number_of_chains = number_of_chains_to_generate;

%hit_record = zeros(numwalkers, 3, number_of_chains);
%capac_record = zeros(numwalkers, 1, number_of_chains);
capacitances = zeros(1,50);
r_infinities = zeros(1,50);

did_it_hit = zeros(50,numwalkers);
chain_time = zeros(1,50);

L=0;

for d = 1:50
    for chain = 1:number_of_chains
    %thickness = 0.002;
    thicknesses(d) = 0.001*d;
    
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
    
    R = sqrt(x_max_mag*x_max_mag + y_max_mag*y_max_mag + z_max_mag*z_max_mag) + 0.5;
    rinfinity = (100 + 20*d)*R;
    r_infinities(d) = rinfinity;
    
    hits = 0;
    losts = 0;

    numsofar = 0;
    
    
    

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
                continue
            elseif r > R
                continue
            else
                for current_plate = 1:number_of_plates_in_chain
                      p = [x y z];
                      %distances(current_plate) = real_distance_1_final(p, collection_of_chains(:,:,current_plate,chain));
                      [in, distances(current_plate)] = new_gjk(p, collection_of_chains(:,:,current_plate,chain));
                      if in == true || (min(distances) < thicknesses(d))
                          hit = true;
                          did_it_hit(d,walkers) = 1;
                          %hit_record(walkers,1,chain) = x;
                          %hit_record(walkers,2,chain) = y;
                          %hit_record(walkers,3,chain) = z;
                          break
                      end
                      %if (min(distances) < 0.02)
                      %    hit = true;
                      %    did_it_hit(chain,walkers) = 1;
                      %    hit_record(walkers,1,chain) = x;
                      %    hit_record(walkers,2,chain) = y;
                      %    hit_record(walkers,3,chain) = z;
                     %     break 
                     % end  
                end
            end
            
            shortstepsize = min(distances);
            
        end
        
    
    end
    
    toc
    chain_time(d) = toc;
    
    for i=1:numwalkers
        capacitances(d) = capacitances(d) + R*did_it_hit(d,i)/numwalkers;
    end
end
end

