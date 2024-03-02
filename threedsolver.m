function u = threedsolver(p,t,f, K, idof, wb, dwb)

% new
w_p = zeros(size(p,2), 1);
if(nargin > 5)
    for i = 1:size(p,2)
        w_p(i) = wb(p(:,i));
    end
else 
    dwb{1} = @(x) 0*x(1)*x(2)*x(3);
    dwb{2} = @(x) 0*x(1)*x(2)*x(3);
    dwb{3} = @(x) 0*x(1)*x(2)*x(3);
end
%

%THREEDSOLVER Summary of this function goes here
%   Detailed explanation goes here

% assumptions: p has three and t has four rows, ie. we use linear elements.
% f is the source function in the Poisson equation.
% K is the geometry matrix.
% idof gives the indices of the interior nodes in the mesh and is of size
% (1,size(p,2))
n_vertices = size(p,2);
Ahat = sparse(n_vertices,n_vertices);
bhat = zeros(n_vertices,1);

% We will use a quadrature rule that gives equal weights to middlepoints of
% the edges on the unit tetrahedron.
load quadrature.mat ti w

% We define the reference basis function on the unit triangle.
phi(:,1) = 1 - ti(1,:) - ti(2,:) - ti(3,:);
phi(:,2) = ti(1,:);
phi(:,3) = ti(2,:);
phi(:,4) = ti(3,:);
dphi = [-1 -1 -1 ;
        1 0 0 ;
        0 1 0 ;
        0 0 1]';

% We loop over the elements, ie. the tetrahedra of the 3D mesh.
for i = 1:size(t,2)
    v = p(:,t(:,i)); % v for vertices of T_i
    sigma = t(:,i);  % sigma gives basis function indices for i:th tetrahedron

    % we define the mapping from the unit tetrahedron to T_i:
    n1 = v(:,1);
    M = [v(:,2)-n1, v(:,3)-n1, v(:,4)-n1];
    b = n1;
    ef = @(x) M*x + b;
    qq = inv(M);

    % We compute contributions to bhat and Ahat:
    for l = 1:4

        for k = 1:4
            Ahat(sigma(l),sigma(k)) = Ahat(sigma(l),sigma(k)) +...
                1/6*dphi(:,l)'*qq*K*qq'*dphi(:,k)*abs(det(M));

            bhat(sigma(l)) = bhat(sigma(l))+...
                abs(det(M))*w(k)*(f(ef(ti(:,k)))*phi(l,k));
            % new
            dwb_Fx = [dwb{1}(ef(ti(:,k))); dwb{2}(ef(ti(:,k))); dwb{3}(ef(ti(:,k)))];
            bhat(sigma(l)) = bhat(sigma(l)) -...
                abs(det(M))*w(k)*(dwb_Fx'*(M'\dphi(:,l)));
            %
        end
    end
end
ndof = n_vertices;
A = Ahat(idof,idof) ; b = bhat(idof);
u = zeros(ndof,1);
u(idof) = A\b;
% new
u = u + w_p;
%