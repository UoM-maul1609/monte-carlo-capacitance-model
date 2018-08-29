function [min_distance] = real_distance_1_final(p,shape_points)

%   a = 1;
%   L = 1;
%   d = 0.5;
%   shape_points = [-a/2 -d*a*sqrt(3)/2 L/2; a/2 -d*a*sqrt(3)/2 L/2; a 0 L/2; a/2 a*sqrt(3)/2 L/2; -a/2 a*sqrt(3)/2 L/2;  -a 0 L/2; -a/2 -d*a*sqrt(3)/2 -L/2; a/2 -d*a*sqrt(3)/2 -L/2; a 0 -L/2; a/2 a*sqrt(3)/2 -L/2; -a/2 a*sqrt(3)/2 -L/2;  -a 0 -L/2];
% % x = 0;
% % y = 0; 
% % z = 0;
% % shape_points = shape_points + [x y z; x y z;x y z;x y z;x y z;x y z;x y z;x y z;x y z;x y z;x y z;x y z];
%  p = [ 1.2417    1.4910    0.4847];

%plot3(position_last(1),position_last(2),position_last(3))
%plot3(p(1),p(2),p(3))
%plot3(shape_points(:,1),shape_points(:,2),shape_points(:,3))

%shape_point = collection_of_chains(:,:,4,10);
simplex = [];
v = p;
farthest = zeros(length(shape_points),1); 

simplexA=[0 0 0];
simplexB=[0 0 0];
simplexC=[0 0 0];

min_distA=0;
min_distB=0;

distance_array = [];
l = 1;

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
                
      
                distance_array(l) = min_distance;
                if (l>1)
                    for k = 1:l-1
                        if(min_distance == distance_array(k))
                            copy_check = 1;
                            min_distance = min(distance_array);
                            break;
                        end
                    end
                end
                l = l+1;
           
                
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
                       
                       perp = cross(cross(simplex(2,:)-simplex(1,:),p-simplex(1,:)),simplex(2,:)-simplex(1,:));
                       v = perp;
                       norm = perp/sqrt(dot(perp,perp));
                       min_distance = dot(p-simplex(1,:),norm);

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
                    
                        perp = cross(CB,CA);                        
                        direction = dot(perp,p-simplex(3,:));

                        if direction > 0                        
                            v = perp;      
                            simplex = [simplex(3,:);simplex(2,:);simplex(1,:)];
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
                           perp = cross(cross(simplex(3,:)-simplex(1,:),p-simplex(1,:)),simplex(3,:)-simplex(1,:));
                           v = perp;
                           norm = perp/sqrt(dot(perp,perp));
                           min_distance = dot(p-simplex(1,:),norm);
                           simplex = [simplex(1,:);simplex(3,:)];
                        else
                            if voro_C_B < 0 & voro_BC > 0 & voro_BC_edge >= 0                                         
                       perp = cross(cross(simplex(3,:)-simplex(2,:),p-simplex(2,:)),simplex(3,:)-simplex(2,:));
                       v = perp;
                       norm = perp/sqrt(dot(perp,perp));
                       min_distance = dot(p-simplex(2,:),norm);
                       simplex = [simplex(2,:);simplex(3,:)];
                            end
                       end
                    end                    
                    end                    
            end
                end
                %fucking tetrahedron
                 if numel(simplex(:,1)) == 4
                    %D_corner
                     D_corner_A = dot(simplex(4,:)-simplex(1,:),p-simplex(4,:));
                     D_corner_B = dot(simplex(4,:)-simplex(2,:),p-simplex(4,:));
                     D_corner_C = dot(simplex(4,:)-simplex(3,:),p-simplex(4,:));
                     
                     if D_corner_A > 0 & D_corner_B > 0 & D_corner_C > 0
                         min_distance = sqrt(dot(simplex(4,:)-p,simplex(4,:)-p));
                         v = p - simplex(4,:);
                         simplex =simplex(4,:);
                         continue
                     end
                     
                     %DC_edge_DCA
                     DC_edge_DCA = cross(cross(simplex(1,:)-simplex(4,:),simplex(3,:)-simplex(4,:)),simplex(3,:)-simplex(4,:));
                     DC_edge_DCA_check = dot(DC_edge_DCA,p-simplex(4,:));
                     
                     CD_edge_BCD = cross(cross(simplex(2,:)-simplex(3,:),simplex(4,:)-simplex(3,:)),simplex(4,:)-simplex(3,:));
                     CD_edge_BCD_check = dot(CD_edge_BCD,p-simplex(4,:));
                     
                     DC_D = dot(simplex(3,:)-simplex(4,:),p-simplex(4,:));
                     DC_C = dot(simplex(4,:)-simplex(3,:),p-simplex(3,:));
                     
                     if DC_edge_DCA_check > 0 & DC_D > 0 & DC_C > 0 & CD_edge_BCD_check > 0
                     perp = cross(cross(simplex(3,:)-simplex(4,:),p-simplex(4,:)),simplex(3,:)-simplex(4,:));
                     norm = perp/sqrt(dot(perp,perp));
                     min_distance = dot(norm,p-simplex(4,:));
                     v = perp;
                     simplex = [simplex(3,:); simplex(4,:)];
                     continue
                     end
                                          
                     
                     %AD_edge_DCA
                     AD_edge_DCA = cross(cross(simplex(3,:)-simplex(1,:),simplex(4,:)-simplex(1,:)),simplex(4,:)-simplex(1,:));
                     AD_edge_DCA_check = dot(AD_edge_DCA,p-simplex(4,:));
                     
                     AD_edge_ABD = cross(cross(simplex(2,:)-simplex(1,:),simplex(4,:)-simplex(1,:)),simplex(4,:)-simplex(1,:));
                     AD_edge_ABD_check = dot(AD_edge_ABD,p-simplex(4,:));
                     
                     AD_A = dot(simplex(4,:)-simplex(1,:),p-simplex(1,:));
                     AD_D = dot(simplex(1,:)-simplex(4,:),p-simplex(4,:));
                     
                     if AD_edge_DCA_check > 0 & AD_D > 0 & AD_A > 0 & AD_edge_ABD_check > 0 
                     perp = cross(cross(simplex(4,:)-simplex(1,:),p-simplex(1,:)),simplex(4,:)-simplex(1,:));
                     norm = perp/sqrt(dot(perp,perp));
                     min_distance = dot(norm,p-simplex(1,:));
                     v = perp;
                     simplex = [simplex(1,:); simplex(4,:)];
                     continue
                     end   
                                              
                     %BD_edge_ABD
                     BD_edge_ABD = cross(cross(simplex(1,:)-simplex(2,:),simplex(4,:)-simplex(2,:)),simplex(4,:)-simplex(2,:));
                     BD_edge_ABD_check = dot(BD_edge_ABD,p-simplex(2,:));
                     
                     BD_edge_BCD = cross(cross(simplex(3,:)-simplex(2,:),simplex(4,:)-simplex(2,:)),simplex(4,:)-simplex(2,:));
                     BD_edge_BCD_check = dot(BD_edge_BCD,p-simplex(2,:));

                     BD_B = dot(simplex(4,:)-simplex(2,:),p-simplex(2,:));
                     BD_D = dot(simplex(2,:)-simplex(4,:),p-simplex(4,:));
                     
                     if BD_edge_ABD_check > 0 & BD_B > 0 & BD_D > 0 & BD_edge_BCD_check > 0
                     perp = cross(cross(simplex(4,:)-simplex(2,:),p-simplex(2,:)),simplex(4,:)-simplex(2,:));
                     norm = perp/sqrt(dot(perp,perp));
                     min_distance = dot(norm,p-simplex(2,:));
                     v = perp;
                     simplex = [simplex(2,:); simplex(4,:)];
                     continue
                     end                                   

                     
                     
                     
%                      %AC_edge_DCA                    
                     AC_edge_DCA = cross(cross(simplex(4,:)-simplex(1,:),simplex(3,:)-simplex(1,:)),simplex(3,:)-simplex(1,:));
                     AC_edge_DCA_check = dot(AC_edge_DCA,p-simplex(1,:));
                     
                     AC_edge_ABC = cross(cross(simplex(2,:)-simplex(3,:),simplex(1,:)-simplex(3,:)),simplex(1,:)-simplex(3,:));
                     AC_edge_ABC_check = dot(AC_edge_ABC,p-simplex(3,:));

                     AC_C = dot(simplex(1,:)-simplex(3,:),p-simplex(3,:));
                     AC_A = dot(simplex(3,:)-simplex(1,:),p-simplex(1,:));
                     
                     if AC_edge_DCA_check > 0 & AC_C > 0 & AC_A > 0 & AC_edge_ABC_check > 0
                     perp = cross(cross(simplex(3,:)-simplex(1,:),p-simplex(3,:)),simplex(3,:)-simplex(1,:));
                     norm = perp/sqrt(dot(perp,perp));
                     min_distance = dot(norm,p-simplex(3,:));
                     v = perp;
                     simplex = [simplex(1,:); simplex(3,:)];
                     continue
                     end
                                                                                                      
                     %DCA face
                     perp = cross(simplex(3,:)-simplex(4,:),simplex(1,:)-simplex(4,:));
                     DCA_face_check = dot(perp,p-simplex(4,:));    
                     if DCA_face_check > 0 & AC_edge_DCA_check <= 0 & AD_edge_DCA_check <= 0 & DC_edge_DCA_check <= 0
                         perp = cross(simplex(3,:)-simplex(4,:),simplex(1,:)-simplex(4,:));
                         norm = perp/sqrt(dot(perp,perp));
                         min_distance = dot(norm,p-simplex(4,:));
                         if min_distance < 0
                             v = -norm;
                         else
                             v = norm;
                         end
                         simplex = [simplex(1,:);simplex(4,:);simplex(3,:)];
                         continue
                     end
                     
                   
                     
%                      %AB_edge_ABD
                     AB_edge_ABD = cross(cross(simplex(4,:)-simplex(1,:),simplex(2,:)-simplex(1,:)),simplex(2,:)-simplex(1,:));
                     AB_edge_ABD_check = dot(AB_edge_ABD,p-simplex(1,:));
                     
                     AB_edge_ABC = cross(cross(simplex(3,:)-simplex(2,:),simplex(1,:)-simplex(2,:)),simplex(1,:)-simplex(2,:));
                     AB_edge_ABC_check = dot(AB_edge_ABC,p-simplex(2,:));
                     
                     AB_A = dot(simplex(2,:)-simplex(1,:),p-simplex(2,:));
                     AB_B = dot(simplex(1,:)-simplex(2,:),p-simplex(1,:));
                     
                     if AB_edge_ABD_check > 0 & AB_A > 0 & AB_B > 0 & AB_edge_ABC_check > 0
                     perp = cross(cross(simplex(2,:)-simplex(1,:),p-simplex(1,:)),simplex(2,:)-simplex(1,:));
                     norm = perp/sqrt(dot(perp,perp));
                     min_distance = dot(norm,p-simplex(1,:));
                     v = perp;
                     simplex = [simplex(1,:); simplex(2,:)];
                     continue
                     end
                                                              
                     %ADB face
                     perp = cross(simplex(1,:)-simplex(4,:),simplex(2,:)-simplex(4,:));
                     ADB_face_check = dot(perp,p-simplex(4,:));
                     if ADB_face_check > 0 & BD_edge_ABD_check <= 0 & AD_edge_ABD_check <= 0 & AB_edge_ABD_check <= 0
                         perp = cross(simplex(4,:)-simplex(1,:),simplex(2,:)-simplex(1,:));
                         norm = perp/sqrt(dot(perp,perp));
                         min_distance = dot(norm,p-simplex(1,:));
                         if min_distance < 0
                             v = -norm;
                         else
                             v = norm;
                         end
                         simplex = [simplex(1,:);simplex(2,:);simplex(4,:)];
                         continue
                     end
                     
                                     
                     

                     
%                      %BC_edge _BCD
                     BC_edge_BCD = cross(cross(simplex(4,:)-simplex(2,:),simplex(3,:)-simplex(2,:)),simplex(3,:)-simplex(2,:));
                     BC_edge_BCD_check = dot(BC_edge_BCD,p-simplex(2,:));
                        
                     BC_edge_ABC = cross(cross(simplex(1,:)-simplex(3,:),simplex(2,:)-simplex(3,:)),simplex(2,:)-simplex(3,:));
                     BC_edge_ABC_check = dot(BC_edge_ABC,p-simplex(3,:));
                     
                     BC_B = dot(simplex(3,:)-simplex(2,:),p-simplex(2,:));
                     BC_C = dot(simplex(2,:)-simplex(3,:),p-simplex(3,:));
                     
                     if BC_edge_BCD_check > 0 & BC_C > 0 & BC_B > 0 & BC_edge_ABC_check > 0
                     perp = cross(cross(simplex(3,:)-simplex(2,:),p-simplex(2,:)),simplex(3,:)-simplex(2,:));
                     norm = perp/sqrt(dot(perp,perp));
                     min_distance = dot(norm,p-simplex(2,:));
                     v = perp;
                     simplex = [simplex(2,:); simplex(3,:)];
                     continue
                     end
                                         
                     %BCD_face
                     perp = cross(simplex(2,:)-simplex(4,:),simplex(3,:)-simplex(4,:));
                     BCD_face_check = dot(perp,p-simplex(4,:));
                     if BCD_face_check > 0 & CD_edge_BCD_check <= 0 & BD_edge_BCD_check <= 0 & BC_edge_BCD_check <= 0
                         perp = cross(simplex(2,:)-simplex(4,:),simplex(3,:)-simplex(4,:));
                         norm = perp/sqrt(dot(perp,perp));
                         min_distance = dot(norm,p-simplex(4,:));
                         if min_distance < 0
                             v = -norm;
                         else
                             v = norm;
                         end
                         simplex = [simplex(4,:);simplex(2,:);simplex(3,:)];
                         continue
                     end                                                                                                                                                                                                                                   
                 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        continue
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

min_distance = sqrt(min_distance^2);
end