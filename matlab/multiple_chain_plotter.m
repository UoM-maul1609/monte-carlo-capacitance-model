rect_face_definitions = [1 2 8 7; 2 3 9 8; 3 4 10 9; 4 5 11 10; 5 6 12 11; 6 1 7 12];
hex_face_definitions = [1 2 3 4 5 6; 7 8 9 10 11 12];

chain_to_plot = 1;

figure('renderer','opengl')
for j=1:number_of_chains_to_generate
    subplot(3,4,j);
    axis equal;

    for i = 1:number_of_plates_in_chain

        patch('Vertices', collection_of_chains(:,:,i,j), 'Faces', rect_face_definitions, 'FaceColor','c')
        patch('Vertices', collection_of_chains(:,:,i,j), 'Faces', hex_face_definitions, 'FaceColor', 'c')

    end
    lighting phong
    light
    view(45,45)
    d=max(D_maxes);
    axis([-d d -d d -d d]./2);
    set(gca,'visible','off');
    axis tight;
%     title({['AI=',num2str(aggregation_indices(j))],...
%         ['Dmax=',num2str(D_maxes(j))],...
%         ['Dmax=',num2str(D_maxes(j))]});
end
