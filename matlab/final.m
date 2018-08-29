%function [min_distance] = real_distance_1(p,shape_points)

 a = 1;
 L = 100;
 shape_points = [-a/2 -a*sqrt(3)/2 L/2; a/2 -a*sqrt(3)/2 L/2; a 0 L/2; a/2 a*sqrt(3)/2 L/2; -a/2 a*sqrt(3)/2 L/2;  -a 0 L/2; -a/2 -a*sqrt(3)/2 -L/2; a/2 -a*sqrt(3)/2 -L/2; a 0 -L/2; a/2 a*sqrt(3)/2 -L/2; -a/2 a*sqrt(3)/2 -L/2;  -a 0 -L/2];
%shape_points = [0.5 0.5 0.5;0.5 0.5 -0.5;0.5 -0.5 -0.5;-0.5 0.5 0.5;-0.5 0.5 -0.5;-0.5 -0.5 -0.5;0.5 -0.5 0.5;-0.5 -0.5 0.5];
p = [3 0 0];

simplex = shape_points(1,:);
min_distance = sqrt(dot(p-simplex(1,:),p-simplex(1,:)));

v = p-simplex;
copy_check = 0;

while copy_check == 0
    
farthest = zeros(length(shape_points),1); 
%support function 1
for i = 1:length(shape_points)
farthest(i) = dot(v,shape_points(i,:))/dot(v,v);
end

[val idx] = (max(farthest));
point = shape_points(idx,:);

last = numel(simplex(:,1))+1;
for i = 1:3
    simplex(last,i) = point(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


                for i = 1:last-1
                if simplex(last,1) == simplex(i,1) && simplex(last,2) == simplex(i,2) && simplex(last,3) == simplex(i,3)
                copy_check = 1;
                break
                else
                copy_check = 0;
                end
                end

                if copy_check == 1
                    break
                end

                if numel(simplex(:,1)) == 2
              
                   AB = simplex(2,:) - simplex(1,:);
                   voro_A = dot(-AB,p-simplex(1,:));
                   voro_B = dot(AB,p-simplex(2,:));
                   
                   if voro_A > 0 
                       v = p-simplex(1,:);     
                       min_distance = sqrt(dot(p-simplex(1,:),p-simplex(1,:)));
                       simplex = simplex(1,:);
                   else                   
                   if voro_B > 0 
                       v = p-simplex(2,:);
                       min_distance = sqrt(dot(p-simplex(2,:),p-simplex(2,:)));
                       simplex = simplex(2,:);
                   else
                       v = -cross(cross(p-simplex(2,:),-AB),-AB);
                       u = dot((p-simplex(2,:)),(simplex(1,:)-simplex(2,:)))/sqrt(dot(simplex(1,:)-simplex(2,:),simplex(1,:)-simplex(2,:)));
                       w = dot((p-simplex(1,:)),(simplex(2,:)-simplex(1,:)))/sqrt(dot(simplex(1,:)-simplex(2,:),simplex(1,:)-simplex(2,:)));
                       tot = u+w;
                       u = u/(tot);
                       w = w/(tot);
                       g = u*simplex(1,:)+w*simplex(2,:);
                       min_distance = sqrt(dot(g-p,g-p));
                   
                   end
                   end
                else
                
            %three points (triangle)    
            if numel(simplex(:,1)) == 3
                    CB = simplex(2,:) - simplex(3,:);
                    CA = simplex(1,:) - simplex(3,:);

                    voro_C_B = dot((-CB),p-simplex(3,:));
                    voro_C_A = dot((-CA),p-simplex(3,:));
                    
                    voro_BC = dot((-CB),p-simplex(2,:));
                    voro_AC = dot((-CA),p-simplex(1,:));
                    
                    AC_edge = cross(cross((simplex(2,:)-simplex(1,:)),-CA),-CA);
                    AB_edge = cross(cross(-CA,simplex(2,:)-simplex(1,:)),simplex(2,:)-simplex(1,:));
                    BC_edge = cross(cross(CA,CB),CB);
                    
                    voro_AC_edge = dot(AC_edge,p-simplex(1,:));
                    voro_AB_edge = dot(AB_edge,p-simplex(2,:));
                    voro_BC_edge = dot(BC_edge,p-simplex(3,:));
                    
                    if voro_AC_edge < 0 & voro_AB_edge < 0 & voro_BC_edge < 0
                    
                        perp = cross(CA,CB);                        
                        direction = dot(perp,p-simplex(3,:));

                        if direction > 0                        
                            v = perp;                            
                        else                            
                            v = -perp;                            
                        end                        
                        norm_perp = perp/sqrt(dot(perp,perp));                        
                        min_distance = dot(p-simplex(3,:),norm_perp); 
                        %simplex = [simplex(2,:);simplex(3,:)];
                  continue
                    else
                    if voro_C_A > 0 & voro_C_B > 0 
                        v = p-simplex(3,:);
                        min_distance = sqrt(dot(simplex(3,:)-p,simplex(3,:)-p));
                        simplex = simplex(3,:);
                    else
                       if voro_C_A < 0 & voro_AC > 0 & voro_AC_edge >= 0
                       v = -cross(cross(p-simplex(3,:),CA),CA);
                       u = dot((p-simplex(3,:)),(simplex(1,:)-simplex(3,:)))/sqrt(dot(simplex(1,:)-simplex(3,:),simplex(1,:)-simplex(3,:)));
                       w = dot((p-simplex(1,:)),(simplex(3,:)-simplex(1,:)))/sqrt(dot(simplex(1,:)-simplex(3,:),simplex(1,:)-simplex(3,:)));
                       tot = u+w;
                       u = u/(tot);
                       w = w/(tot);
                       g = u*simplex(1,:)+w*simplex(3,:);
                       min_distance = sqrt(dot(g-p,g-p));
                       simplex = [simplex(1,:);simplex(2,:)];
                        else
                            if voro_C_B < 0 & voro_BC > 0 & voro_BC_edge >= 0 
                       v = -cross(cross(p-simplex(3,:),CB),CB);
                       u = dot((p-simplex(3,:)),(simplex(2,:)-simplex(3,:)))/sqrt(dot(simplex(2,:)-simplex(3,:),simplex(2,:)-simplex(3,:)));
                       w = dot((p-simplex(2,:)),(simplex(3,:)-simplex(2,:)))/sqrt(dot(simplex(2,:)-simplex(3,:),simplex(2,:)-simplex(3,:)));
                       tot = u+w;
                       u = u/(tot);
                       w = w/(tot);
                       g = u*simplex(2,:)+w*simplex(3,:);
                       min_distance = sqrt(dot(g-p,g-p));
                       simplex = [simplex(2,:);simplex(3,:)];
                            end
                       end
                    end                    
                    end                    
            end
                end
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        continue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

min_distance = sqrt(min_distance^2)
%end