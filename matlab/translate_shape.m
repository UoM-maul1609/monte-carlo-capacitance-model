function [ translated_shape_points ] = translate_shape( shape_points, dx, dy, dz )

translated_shape_points(:,1) = shape_points(:,1) + dx;
translated_shape_points(:,2) = shape_points(:,2) + dy;
translated_shape_points(:,3) = shape_points(:,3) + dz;

end

