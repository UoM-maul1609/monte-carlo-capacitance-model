R = 4;
rinfinity = 500*R;
numwalkers = 10000;

number_of_chains = number_of_chains_to_generate;

hit_record = zeros(numwalkers, 3, number_of_chains);
capac_record = zeros(numwalkers, 1, number_of_chains);
capacitances = zeros(1,number_of_chains);



for chain = 1:number_of_chains

    shortstepsize = 0.01;
    hits = 0;
    losts = 0;

    numsofar = 0;
    
    x_max = max(max(collection_of_chains(:,1,:,chain)));
    x_min = min(min(collection_of_chains(:,1,:,chain)));
    
    y_max = max(max(collection_of_chains(:,2,:,chain)));
    y_min = min(min(collection_of_chains(:,2,:,chain)));
    
    z_max = max(max(collection_of_chains(:,3,:,chain)));
    z_min = min(min(collection_of_chains(:,3,:,chain)));

    parfor walkers = 1:numwalkers
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
                losts = losts + 1;
                continue
            elseif r > R
                continue
            else
                if(x > x_max || x < x_min || y > y_max || y < y_min || z > z_max || z < z_min)
                    continue;
                end
                for current_plate = 1:number_of_plates_in_chain
                      p = [x y z];
                      in = new_gjk(p, collection_of_chains(:,:,current_plate,chain));
                      if(in == true)
                          hit = true;
                          hits = hits+1;
                          hit_record(walkers,1,chain) = x;
                          hit_record(walkers,2,chain) = y;
                          hit_record(walkers,3,chain) = z;
                          break; 
                      end  
                end
            end
        end

        chain
        walkers
        ratio = hits/(hits+losts)
        capacitance = (hits / (hits + losts))*R
        capac_record(walkers,chain) = (hits / (hits + losts))*R;
    
    end

    capacitances(chain) = (hits / (hits + losts))*R;
end

