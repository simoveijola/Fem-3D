function outer_nodes = solve_outer_nodes(t)

    n = size(t,2);
    faces = zeros(3, 4*n);

    faceCounts = containers.Map('KeyType','int64','ValueType','int64');
    keys = zeros(3,size(faces,2));

    for i = 1:n
        
        nodes = [t(1,i); t(2,i); t(3,i); t(4,i)];

        faces(:, 4*i + 1) = [nodes(1); nodes(2); nodes(3)];
        faces(:, 4*i + 2) = [nodes(1); nodes(2); nodes(4)];
        faces(:, 4*i + 3) = [nodes(1); nodes(3); nodes(4)];
        faces(:, 4*i + 4) = [nodes(2); nodes(3); nodes(4)];

        for j = 1:4
            face = sort(faces(:,4*i + j));
            key = hashCode(face);
            keys(:, 4*i + j) = face;
            if(isKey(faceCounts, key)) 
                faceCounts(key) = faceCounts(key) + 1;
            else
                faceCounts(key) = 1;
            end
        end

    end

    keys = setdiff(unique(keys','rows'), [0,0,0], 'rows')';
    outer_nodes = [];

    for i = 1:length(keys)
        key = hashCode(keys(:,i));
        if(faceCounts(key) == 1)
            outer_nodes = union(outer_nodes, keys(:,i));
        end
    end
    
end

%% 
function h = hashCode(a)

    if(length(a) ~= 3)
        disp('asd')
    end
    P1 = 139;
    P2 = 163;

    h = (a(1)*P1 + a(2))*P2 + a(3);

end