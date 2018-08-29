function [in] = new_gjk(p, shape_points)


v = p;

simplex = [];
inside = -1;
farthest = zeros(length(shape_points),1); 

%support function
    for i = 1:length(shape_points)
          dot = v(1)*shape_points(i,1) + v(2)*shape_points(i,2) + v(3)*shape_points(i,3);
        farthest(i) = dot;
    end

    [val idx] = (max(farthest));
    point_1 = shape_points(idx,:);
    point_3 = point_1 - p;

        for i = 1:3
            simplex(i) = point_3(i);
        end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first check to see if P is outside the shape
    check = point_3(1)*p(1) + point_3(2)*p(2) + point_3(3)*p(3);
    if check < 0;
        inside = 0;
        in = false;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = -v;

%loop for adding new points
while(inside == -1)

%second support function
    dot = 0;
    for i = 1:length(shape_points)
        for j = 1:3
        dot = dot + v(j)*shape_points(i,j);
        end
        farthest(i) = dot;
        dot = 0;
    end

    [val idx] = (max(farthest));
    point_1 = shape_points(idx,:);
    point_3 = point_1 - p;

        last = numel(simplex(:,1))+1;
        for i = 1:3
        simplex(last,i) = point_3(i);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%check the last simplex point is within the direction of the origin
    last = numel(simplex(:,1));
    
    dot = 0;
    for i = 1:3
    dot = simplex(last,1)*v(1) + simplex(last,2)*v(2) + simplex(last,3)*v(3);
    %dot = dot + simplex(last,i)*v(i); 
    end
    
    if dot < 0;
        inside = 0;
        in = false;
    continue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else       
            %check for origin   
            %two points (line)
                if numel(simplex(:,1)) == 2
                    AB = simplex(2,:) - simplex(1,:);

                    simp1_AB = 0;
                    AB_AB = 0;
                    for i = 1:3
                        simp1_AB = simp1_AB + simplex(1,i)*AB(i);
                        AB_AB = AB_AB + AB(i)*AB(i);
                    end      

                    
                   v = AB*(simp1_AB) - simplex(1,:)*(AB_AB);
                      if v(1) == 0 & v(2) == 0 & v(3) == 0
                          inside = 1;
                      end
                else
                
            %three points (triangle)    
                      if numel(simplex(:,1)) == 3
                    CB = simplex(2,:) - simplex(3,:);
                    CA = simplex(1,:) - simplex(3,:);

                    CB_CA = 0;
                    CB_CB = 0;
                    CA_CA = 0;
                    for i = 1:3
                        CB_CA = CB_CA + CA(i)*CB(i);
                        CB_CB = CB_CB + CB(i)*CB(i);
                        CA_CA = CA_CA + CA(i)*CA(i);
                    end

                    region_1 = CB*(CB_CA) - CA*(CB_CB);
                    region_2 = CA*(CB_CA) - CB*(CA_CA);

                    origin_1 = 0;
                    origin_2 = 0;
                    for i = 1:3
                        origin_1 = origin_1 - simplex(3,i)*region_1(i);
                        origin_2 = origin_2 - simplex(3,i)*region_2(i);
                    end
                    
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
                             perp = zeros(1,3);
                             perp(1) = CB(2)*CA(3) - CB(3)*CA(2);
                             perp(2) = CB(3)*CA(1) - CB(1)*CA(3);
                             perp(3) = CB(1)*CA(2) - CB(2)*CA(1);

                             origin_3 = 0;
                             for i = 1:3
                                 origin_3 = origin_3 - perp(i)*simplex(3,i);
                             continue
                             end

                              if origin_3 == 0
                                  inside = 1;
                              end
                             
                             if origin_3 > 0
                                 v = perp;
                                 simp = [simplex(3,:);simplex(2,:);simplex(1,:)];
                                 simplex = simp;
                             continue
                             end
                             
                             if origin_3 < 0
                                 v = -perp;
                                 %simplex = [simplex(2,:);simplex(1,:);simplex(3,:)];
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
                
                                                
                        perp_DAB = zeros(1,3);
                        perp_DBC = zeros(1,3);
                        perp_DCA = zeros(1,3);

                    perp_DAB(1,1) = DA(2)*DB(3) - DA(3)*DB(2);
                    perp_DAB(1,2) = DA(3)*DB(1) - DA(1)*DB(3);
                    perp_DAB(1,3) = DA(1)*DB(2) - DA(2)*DB(1);

                    perp_DBC(1,1) = DB(2)*DC(3) - DB(3)*DC(2);
                    perp_DBC(1,2) = DB(3)*DC(1) - DB(1)*DC(3);
                    perp_DBC(1,3) = DB(1)*DC(2) - DB(2)*DC(1);

                    perp_DCA(1,1) = DC(2)*DA(3) - DC(3)*DA(2);
                    perp_DCA(1,2) = DC(3)*DA(1) - DC(1)*DA(3);
                    perp_DCA(1,3) = DC(1)*DA(2) - DC(2)*DA(1);
                  
                    origin_DAB = 0;
                    origin_DBC = 0;
                    origin_DCA = 0;

                    for i = 1:3
                        origin_DAB = origin_DAB - perp_DAB(i)*simplex(4,i);
                        origin_DBC = origin_DBC - perp_DBC(i)*simplex(4,i);
                        origin_DCA = origin_DCA - perp_DCA(i)*simplex(4,i);
                    end

                    %origin is within the tetrahedron
                    if origin_DAB <= 0 & origin_DBC <= 0 & origin_DCA <= 0
                        inside = 1;
                    end

                   
                    
                    if origin_DAB > 0
                        simp = [simplex(1,:); simplex(2,:); simplex(4,:)];
                        simplex = simp;
                        v = perp_DAB;
                        continue
                    end
                    
                    if origin_DBC > 0
                        simp = [simplex(4,:);simplex(2,:); simplex(3,:)];
                        simplex = simp;
                        v = perp_DBC;
                        continue
                    end
                    
                    if origin_DCA > 0
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
        

if inside == 1
    in = true;
end

%in
end