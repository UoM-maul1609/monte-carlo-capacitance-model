function simplex = support_function(v,p,shape_points,simplex)


farthest = zeros(length(shape_points),1);  
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
point_2 = p;
point_3 = point_1 - point_2;

if numel(simplex)==0
    for i = 1:3
        simplex(i) = point_3(i);
    end
else
    last = numel(simplex(:,1))+1;
    for i = 1:3
    simplex(last,i) = point_3(i);
    end
end
end