function in = new_gjk_polyhedra(shape_points_1, shape_points_2)

v = [4 0 0];
simplex = [];
inside = -1;
farthest_1 = zeros(size(shape_points_1,1),1); 
farthest_2 = zeros(size(shape_points_2,1),1);
%support function
    for i = 1:size(shape_points_1,1)
        farthest_1(i) = dot(v,shape_points_1(i,:));
    end

    for i = 1:size(shape_points_2,1)       
        farthest_2(i) = dot(-v,shape_points_2(i,:));
    end
    
    [val_1 idx_1] = (max(farthest_1));
    [val_2 idx_2] = (max(farthest_2));
    point_1 = shape_points_1(idx_1,:);
    point_2 = shape_points_2(idx_2,:);
    point_3 = point_1 - point_2;

        for i = 1:3
            simplex(i) = point_3(i);
        end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = -v;

%loop for adding new points
while(inside == -1)

%second support function
    for i = 1:size(shape_points_1,1)
        farthest_1(i) = dot(v,shape_points_1(i,:));
    end

    for i = 1:size(shape_points_2,1)
        farthest_2(i) = dot(-v,shape_points_2(i,:));
    end
    
    [val_1 idx_1] = (max(farthest_1));
    [val_2 idx_2] = (max(farthest_2));
    point_1 = shape_points_1(idx_1,:);
    point_2 = shape_points_2(idx_2,:);
    point_3 = point_1 - point_2;
    
        last = numel(simplex(:,1))+1;
        for i = 1:3
        simplex(last,i) = point_3(i);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%check the last simplex point is within the direction of the origin
    last = numel(simplex(:,1));
    check = dot(simplex(last,:),v);

    if check < 0;
        inside = 0;
        in = false;
    continue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else       
            %check for origin   
            %two points (line)
                if numel(simplex(:,1)) == 2
                    AB = simplex(2,:) - simplex(1,:);
                    
                   v = cross(AB,cross(AB,simplex(1,:)));
                      continue
                else
                
            %three points (triangle)    
                      if numel(simplex(:,1)) == 3
                    CB = simplex(2,:) - simplex(3,:);
                    CA = simplex(1,:) - simplex(3,:);
                   
                    region_1 = cross(cross(CA,CB),CB);
                    region_2 = cross(cross(CB,CA),CA);

                    origin_1 = dot(-simplex(3,:),region_1);
                    origin_2 = dot(-simplex(3,:),region_2);
                    
                    if origin_1 > 0
                         v = region_1;    
                         %get rid of A
                         simplex = [simplex(2,:);simplex(3,:)];
                    continue 
                    end

                     if origin_2 > 0
                         v = region_2;
                         %remove point B
                         simplex = [simplex(1,:);simplex(3,:)];
                     continue
                     end
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %origin is within triangle, now find direction to make tetrahedron
                       if origin_1 && origin_2 < 0

                             %perp vector to the triangle
                             perp = cross(CB,CA);

                             origin_3 = dot(perp,simplex(3,:));

                             if origin_3 == 0
                                  inside = 1;
                                  continue
                             end
                             
                             if origin_3 > 0
                                 v = perp;
                                 simp = [simplex(3,:);simplex(2,:);simplex(1,:)];
                                 simplex = simp;
                             continue
                             end
                             
                             if origin_3 < 0
                                 v = -perp;
                                continue
                             end
                        end
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      else

                    %four points (tetrahedron)
                    if numel(simplex(:,1)) == 4
                        DC = simplex(3,:) - simplex(4,:);
                        DB = simplex(2,:) - simplex(4,:);
                        DA = simplex(1,:) - simplex(4,:);
                                   
                    perp_DAB = cross(DA,DB);

                    perp_DBC = cross(DB,DC);

                    perp_DCA = cross(DC,DA);

                  
                    origin_DAB = dot(-perp_DAB,simplex(4,:));
                    origin_DBC = dot(-perp_DBC,simplex(4,:));
                    origin_DCA = dot(-perp_DCA,simplex(4,:));

                    %origin is within the tetrahedron
                    if origin_DAB <= 0 & origin_DBC <= 0 & origin_DCA <= 0
                        inside = 1;
                        continue
                    end
                   
                    if origin_DAB >= 0
                        simp = [simplex(1,:); simplex(2,:); simplex(4,:)];
                        simplex = simp;
                        v = perp_DAB;
                        continue
                    end
                    
                    if origin_DBC >= 0
                        simp = [simplex(4,:);simplex(2,:); simplex(3,:)];
                        simplex = simp;
                        v = perp_DBC;
                        continue
                    end
                    
                    if origin_DCA >= 0
                        simp = [simplex(1,:); simplex(4,:); simplex(3,:)];
                        simplex = simp;
                        v = perp_DCA;
                        continue
                    end
                    end
                      end
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end of origin check

    end
%end of dot test

end
%end of while loop    
if inside == 0 
in = false;
end

if inside == 1
    in = true;
end
inside
end