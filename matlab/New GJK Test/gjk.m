function [in] = gjk(p, shape_points)

v = p;
simplex = [];
inside = -1;
%support function
simplex = support_function(v,p,shape_points,simplex);
v = -v;

while(inside == -1)
%point b, second element of the simplex
simplex = support_function(v,p,shape_points,simplex);

%check the last simplex point is within the direction of the origin
last = numel(simplex(:,1));
dot = 0;
for i = 1:3
    dot = dot + simplex(last,i)*v(i);
end
    if dot < 0;
        inside = 0;
        in = false;
    continue
    else    
       
%check for origin   
%inside = origin_check(simplex);

[v, simplex, inside] = origin_check(v, simplex, inside);


continue
    end
    
end     

if inside == 1
    in = true;
end


end
    
