
function [x_out, t_out, U_out] = diffusion1d( kappa, x_rng, n_x, t_rng, n_t, u_init, u_bndry )
%-----------------------------------------------
% Diffuion Equation by Mihir Kakkar
%-----------------------------------------------
% Parameters
% ==========
%    kappa This is the diffusivity coefficient
%    x_rng Array describing the spatial boundaries
%    t_rng Array describing the initial and final times
%
%    u_init This describes the initial distribution of heat across space at t=0
%    u_bndry This describes the u-value at each of the boundaries at a given time
%
%    n_x Number of points along the x axis, incrementing from a to b
%    n_t Number of points along the t axis, incrementing from t_initial, t_final
%-----------------------------------------------
% Return Values
% =============
%    x_out The vector x is an array of n_x x values , going from a to b
%    t_out The vector t is an array of n_t x values , going from t_initial, t_final
%    U_out The matrix U describes the U values at each point in the n_x by n_t matrix
%-----------------------------------------------

% Error checking
 
    if ~all( size( x_rng ) == [2, 1] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument x_rng is not a 2-dimensional column vector' ) );
    end
 
    if ~all( size( t_rng  ) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument t_rng is not a 2-dimensional column vector' ) );
    end
  
    if ~isscalar( kappa ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument kappa is not a scalar' ) );
    end
 
    if ~isscalar( n_x ) || ( n_x ~= round( n_x ) )  
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument n_x is not an integer' ) );
    end
 
    if ~isscalar( n_t ) || ( n_t ~= round( n_t ) )  
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument n_t is not an integer' ) );
    end
 
    if ~isa( u_init, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_init is not a function handle' ) );
    end
 
    if ~isa( u_bndry, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_bndry is not a function handle' ) );
    end
    
    h = (x_rng(2)- x_rng(1))/(n_x - 1);
    dt = (t_rng(2)-t_rng(1))/(n_t - 1);
    
    error_check = sprintf('The ratio kappa * dt / h^2 = %d, consider using n_t = %d', kappa * dt / h^2, floor(kappa * (t_rng(2)-t_rng(1))/(0.5 * h^2) + 2) );
    
    if kappa * dt / h^2 >= 0.5
        throw ( MException('MATLAB:invalid_argument', ... 
        error_check));
    end
    
% Initialization of x and t vectors, and matrix M
 
    x_out = linspace(x_rng(1), x_rng(2), n_x)';
    t_out = linspace(t_rng(1), t_rng(2), n_t);
 
    M = zeros(n_x, n_t);
    u = u_bndry(t_out);
 
    M(1:end,1) = u_init(x_out);
    M(1, 2:end) = u(1, 2:end);
    M(end, 2:end) = u(2, 2:end);
 
% Solve
 
    for i = 1:n_t - 1 
        u_i = M(:,i);
        M(2:end-1, i+1) = u_i(2:end-1) +  diff(u_i, 2) * kappa * dt / h^2;
    end
    
    U_out = M(:,:);
 
end


