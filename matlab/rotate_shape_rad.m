function [ rotated_points ] = rotate_shape_rad( shape_points, theta_x_rad, theta_y_rad, theta_z_rad )
% This function rotates the points of shape_points first about x, then
% about y, then about z by the angles theta_x, theta_y, theta_z.


rotated_points(:,1) = (cos(theta_y_rad)*cos(theta_z_rad))*shape_points(:,1) + (cos(theta_z_rad)*sin(theta_x_rad)*sin(theta_y_rad) - cos(theta_x_rad)*sin(theta_z_rad))*shape_points(:,2) + (cos(theta_x_rad)*cos(theta_z_rad)*sin(theta_y_rad) + sin(theta_x_rad)*sin(theta_z_rad))*shape_points(:,3);
rotated_points(:,2) = (cos(theta_y_rad)*sin(theta_z_rad))*shape_points(:,1) + (cos(theta_x_rad)*cos(theta_z_rad) + sin(theta_x_rad)*sin(theta_y_rad)*sin(theta_z_rad))*shape_points(:,2) + (-cos(theta_z_rad)*sin(theta_x_rad) + cos(theta_x_rad)*sin(theta_y_rad)*sin(theta_z_rad))*shape_points(:,3);
rotated_points(:,3) = (-sin(theta_y_rad))*shape_points(:,1) + (cos(theta_y_rad)*sin(theta_x_rad))*shape_points(:,2) + (cos(theta_x_rad)*cos(theta_y_rad))*shape_points(:,3);

end