function [U_out] = laplace2d( U )
    [n_x, n_y] = size( U );
    U_out = U;
 
% Error and Warning Checking
% ==========================
%
% Ensuring boundaries are all defined at edges, and that the inputted matrix is 2D
% U must be a 2D matrix
if ndims(U) ~= 2
throw( MException( 'MATLAB:invalid_argument', 'this is not a 2D  matrix' ) )
end

% Boundary conditions are defined at the edges
inf_check_vert = isinf(-U(1, :)) + isinf(-U(end, :));
inf_check_horiz = isinf(-U(:, 1)) + isinf(-U(:, end));
    
if max(inf_check_vert) >= 1 || max(inf_check_horiz) >= 1
   throw( MException( 'MATLAB:invalid_argument', 'Boudary Conditions must be defined at the Boundaries' ) ) 
    end
% Initialization
% ==============
%
% Create an empty solution matrix U_out with same dimensions as the input matrix. %Create u_to_w and w_to_u matrices.
  
  M = spalloc( m, m, 5*m );
    b = zeros( m, 1 );

% Mapping the unknown points to a unique number from 1 to m
% =========================================================
%
% Create a system of equations over all the points to find unknown values.


    for ix = 1:n_x
        for iy = 1:n_y
            if U(ix, iy) == -Inf
                m = m + 1;
                u_to_w(ix, iy) = m;
                w_to_u(:, m) = [ix, iy]';
            end
        end
    end

% Creating and solving a system of linear equations
% =================================================
%
%   Your description here.

    for k = 1:m
        % Get the coordinates of the kth point
        k_coords = w_to_u(:,k);
        % For each of the 4 adjacent points, determine if
        % the point is an insluated boundary point, a Dirichlet
        % boundary point or an unknown value and modify M as appropriate.
        for y_offset = [-1, 1]
            nbr_coords = [k_coords(1) + y_offset, k_coords(2)];
            nbr_U = U(nbr_coords(1),nbr_coords(2)); 
            if isnan(nbr_U)
                continue
            elseif nbr_U == -Inf
                M(k,k) = M(k,k) - 1;               
                nbr_k = u_to_w(nbr_coords(1),nbr_coords(2));
                M(k, nbr_k) = M(k, nbr_k) + 1; 
            else 
                M(k,k) = M(k,k) - 1;
                b(k) = b(k) - nbr_U;
            end
        end
        for x_offset = [-1, 1]
            nbr_coords = [k_coords(1), k_coords(2) + x_offset];
            nbr_U = U(nbr_coords(1),nbr_coords(2)); 
            if isnan(nbr_U)
                continue
            elseif nbr_U == -Inf
                M(k,k) = M(k,k) - 1;               
                nbr_k = u_to_w(nbr_coords(1),nbr_coords(2));
                M(k, nbr_k) = M(k, nbr_k) + 1; 
            else 
                M(k,k) = M(k,k) - 1;
                b(k) = b(k) - nbr_U; 
            end
        end
    end

% Substituting the values back into the matrix U_out
% ===================================================
%
%   w_to_u provides the location of each point in the matrix being iterated. Using this, %     each of the calculated values can be put into that coordinate in the U_out solution %matrix.

    w = M \ b;
 
    for k = 1:m
        k_coord = w_to_u(:,k);
        U_out(k_coord(1), k_coord(2)) = w(k);
    end
end




 

 

