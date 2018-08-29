R = 1;
rinfinity = 500*R;
numwalkers = 1e5;

shape_points = [0.5 0.5 0.5;0.5 0.5 -0.5;0.5 -0.5 -0.5;-0.5 0.5 0.5;-0.5 0.5 -0.5;-0.5 -0.5 -0.5;0.5 -0.5 0.5;-0.5 -0.5 0.5];
epsilon = 0.001;

shortstepsize = 0.01;
hits = 0;
losts = 0;

numsofar = 0;

%for hit record
hit_record = zeros(numwalkers,3);
%

%for capacitance_record
capac_record = zeros(numwalkers,1);
%


for walkers = 1:numwalkers
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
            p = [x y z];
            in = new_gjk(p, shape_points);
            %in = gjk(p,shape_points);
           if in == true
                hit = true;
                hits = hits + 1;
             
                %record the position of the hits to make sure gjk isn't
                %fucked
                hit_record(walkers,1) = x;
                hit_record(walkers,2) = y;
                hit_record(walkers,3) = z;
                %
                
                continue
            else
                continue
            end
        
        end
    end
    
    walkers
    ratio = (hits / (hits + losts))
    capacitance = (hits / (hits + losts))*R
    
    capac_record(walkers) = capacitance;
    
end
        
        
        
        
        