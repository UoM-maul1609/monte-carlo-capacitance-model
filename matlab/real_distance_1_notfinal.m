function [min_distance] = real_distance_1(p,shape_points)

% p = [0.3   0.1   -2];
% a = 1;
% L = 1;
% shape_points = [-a/2 -a*sqrt(3)/2 L/2; a/2 -a*sqrt(3)/2 L/2; a 0 L/2; a/2 a*sqrt(3)/2 L/2; -a/2 a*sqrt(3)/2 L/2;  -a 0 L/2; -a/2 -a*sqrt(3)/2 -L/2; a/2 -a*sqrt(3)/2 -L/2; a 0 -L/2; a/2 a*sqrt(3)/2 -L/2; -a/2 a*sqrt(3)/2 -L/2;  -a 0 -L/2];
%shape_points = [0.5 0.5 0.5;0.5 0.5 -0.5;0.5 -0.5 -0.5;-0.5 0.5 0.5;-0.5 0.5 -0.5;-0.5 -0.5 -0.5;0.5 -0.5 0.5;-0.5 -0.5 0.5];


simplex = [];
v = p;
farthest = zeros(length(shape_points),1); 
%support function 1
for i = 1:length(shape_points)
farthest(i) = dot(v,shape_points(i,:))/dot(v,v);
end

[val idx] = (max(farthest));
point = shape_points(idx,:);

for i = 1:3
    simplex(i) = point(i);
end

%simplex = shape_points(1,:);
min_distance = sqrt(dot(p-simplex(1,:),p-simplex(1,:)));
v = p-simplex;
copy_check = 0;

%begin loop
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
                       simplex = [simplex(1,:);simplex(3,:)];
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
                
%                  if numel(simplex(:,1)) == 4
%                     %BDA    
%                     perp_1 = cross(simplex(4,:)-simplex(2,:),simplex(1,:)-simplex(2,:));                        
%                         direction_1 = dot(perp_1,p-simplex(2,:));
% 
%                         if direction_1 > 0                        
%                             v_1 = perp_1;                            
%                         else                            
%                             v_1 = -perp_1;                            
%                         end                        
%                         norm_perp_1 = perp_1/sqrt(dot(perp_1,perp_1));                        
%                         min_distance_1 = dot(p-simplex(2,:),norm_perp_1);              
%                     
%                     %BDC   
%                     perp_2 = cross(simplex(2,:)-simplex(4,:),simplex(3,:)-simplex(4,:));                        
%                         direction_2 = dot(perp_2,p-simplex(4,:));
% 
%                         if direction_2 > 0                        
%                             v_2 = perp_2;                            
%                         else                            
%                             v_2 = -perp_2;                            
%                         end                        
%                         norm_perp_2 = perp_2/sqrt(dot(perp_2,perp_2));                        
%                         min_distance_2 = dot(p-simplex(4,:),norm_perp_2);              
%                     
%                      %ADC
%                      perp_3 = cross(simplex(4,:)-simplex(1,:),simplex(3,:)-simplex(1,:));                        
%                         direction_3 = dot(perp_3,p-simplex(1,:));
% 
%                         if direction_3 > 0                        
%                             v_3 = perp_3;                            
%                         else                            
%                             v_3 = -perp_3;                            
%                         end                        
%                         norm_perp_3 = perp_3/sqrt(dot(perp_3,perp_3));                        
%                         min_distance_3 = dot(p-simplex(1,:),norm_perp_3); 
%                     
%                         
%                      [value_min idx_min] = min([min_distance_1 min_distance_2 min_distance_3]);
%                      min_distance = value_min;
%                      if idx_min == 1
%                          v = v_1;
%                          simplex = [simplex(1,:);simplex(2,:);simplex(4,:)];
%                      end
%                      
%                      if idx_min == 2
%                          v = v_2;
%                          simplex = [simplex(2,:);simplex(3,:);simplex(4,:)];
%                      end
%                      
%                      if idx_min == 3
%                          v = v_3;
%                          simplex = [simplex(1,:);simplex(3,:);simplex(4,:)];
%                      end
%                 end
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        continue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

min_distance = sqrt(min_distance^2);
end