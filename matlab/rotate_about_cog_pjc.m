% rotate about COG
plot1=false;
if k==1
    clear D_m1 area1 v1;
end
rect_face_definitions = [1 2 8 7; 2 3 9 8; 3 4 10 9; 4 5 11 10; 5 6 12 11; 6 1 7 12];
hex_face_definitions = [1 2 3 4 5 6; 7 8 9 10 11 12];
if plot1
    figure('renderer','openGl');;
end
chain_to_plot = k;
theta1=rand(1,1).*2.*pi;
theta2=rand(1,1).*2.*pi;
theta3=rand(1,1).*2.*pi;
Rx=[1 0 0;0 cos(theta1) -sin(theta1); 0 sin(theta1) cos(theta1)];
Ry=[cos(theta2) 0 sin(theta2); 0 1 0;-sin(theta2) 0 cos(theta2)];
Rz=[cos(theta3) -sin(theta3) 0; sin(theta3) cos(theta3) 0; 0 0 1];
R=Rx*Ry*Rz;
for i=1:number_of_plates_in_chain
    for j=1:12
        pos(j,:,i)=R*squeeze(collection_of_chains(j,:,i,chain_to_plot))';
    end
end

if plot1
    for i = 1:number_of_plates_in_chain

        h=patch('Vertices', pos(:,1:2,i,chain_to_plot), 'Faces', rect_face_definitions, 'FaceColor','c')
        set(h,'facealpha',0.5);
        h=patch('Vertices', pos(:,1:2,i,chain_to_plot), 'Faces', hex_face_definitions, 'FaceColor', 'c')
        set(h,'facealpha',0.5);

    end
    axis equal;
end


% find the union of all the polygons;
xa=pos(hex_face_definitions(1,:),1,1);
ya=pos(hex_face_definitions(1,:),2,1);
for i=1:number_of_plates_in_chain
    for j=1:6
        [xa,ya]=polybool('union',xa,...
            ya,...
            pos(rect_face_definitions(j,:),1,i),...
            pos(rect_face_definitions(j,:),2,i) );
    end
    for j=1:2
        [xa,ya]=polybool('union',xa,...
            ya,...
            pos(hex_face_definitions(j,:),1,i),...
            pos(hex_face_definitions(j,:),2,i) );
    end
end

ind=find(isnan(xa));
while length(ind)>0
    [xa,ya]=polybool('union',...
        xa(1:ind(1)-1),...
        ya(1:ind(1)-1),...
        xa(ind(1)+1:end),...
        ya(ind(1)+1:end) );
     ind=find(isnan(xa));   
end

if plot1
    hold on;
    plot(xa,ya,'m');
end

clear all1;
area=polyarea(xa,ya);
for i=1:length(xa)
    for j=1:length(xa)
       all1(i,j)=sqrt((xa(i)-xa(j)).^2+(ya(i)-ya(j)).^2); 
    end
end
D_m=max(all1(:));

unit1=40e-6;
area1(k)=area.*unit1.^2;
if isnan(area)
    area
end
D_m1(k)=D_m.*unit1;
v1(k)=number_of_plates_in_chain.*3.*sqrt(3)./2.*a.^2.*0.2.*unit1.^3;
Ar=area1./(pi./4.*D_m1.^2);
eta=2e-5;

rhoa=20000./287./210.;
X=rhoa.*8.*v1.*910.*9.8./(eta.^2.*pi.*Ar.^0.5);
Re=8.^2./4.*((1+4.*sqrt(X)./(8.^2.*sqrt(0.35))).^0.5-1).^2;
vt=eta.*Re./(rhoa.*D_m1);

