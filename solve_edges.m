function outer_nodes = solve_edges(p,t)

    n = size(t,2);
    faces = zeros(3, 4*n);
    
    faceCounter = zeros(size(faces,2),1);

    for i = 1:n
        
        nodes = [t(1,i); t(2,i); t(3,i); t(4,i)];

        faces(:, 4*i + 1) = [nodes(1); nodes(2); nodes(3)];
        faces(:, 4*i + 2) = [nodes(1); nodes(2); nodes(4)];
        faces(:, 4*i + 3) = [nodes(1); nodes(3); nodes(4)];
        faces(:, 4*i + 4) = [nodes(2); nodes(3); nodes(4)];
        faceCount 

    end
    
    for i = 1:size(faces,2)
        for j = 1:size(faces,2)
            if(length(intersection(faces(:,i), faces(:,j))) == 3)
                faceCounter(i) = faceCounter(i) + 1;
            end
        end
    end

    I = find(faceCounter == 1);
    outer_faces = faces(:,I);
    outer_nodes = unique(outer_faces);

end