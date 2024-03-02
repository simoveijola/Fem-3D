function p = plot3D_solution(p,t,uh)

    X = zeros(4, size(t,2));
    Y = zeros(4, size(t,2));
    Z = zeros(4, size(t,2));
    C = zeros(4, size(t,2));

    for i = 1:size(t,2)
        X(:,i) = p(1, t(:,i));
        Y(:,i) = p(2, t(:,i));
        Z(:,i) = p(3, t(:,i));
        C(:,i) = uh(t(:,i));
    end

    p = patch(X,Y,Z,C,'faceAlpha', 0.05, 'LineStyle', 'none');
    colormap(jet(256))
    colorbar
    
end