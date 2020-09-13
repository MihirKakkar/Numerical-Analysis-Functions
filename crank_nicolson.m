function [x_out, t_out, U_out] = crank_nicolson( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )

%-----------------------------------------------
% Crank Nicolson by Mihir Kakkar
%-----------------------------------------------
% Crank Nicolson 1D with Neumann Conditoin Apllied to boundary B: 
% An unconditionally stable method of solvning 1D equations
%-----------------------------------------------
% Inputs: 
%-----------------------------------------------
%   kappa - Coefficient of Diffusion
%	x_rng - Array of range over which space is defined
%	nx - Step size of space
%	t_rng - Array of range over which time is defined
%	nt - Step size of time
%	u_init - function handle for solving inital values
%	u_bndry - finction handle for solving boundary values
%-----------------------------------------------
% Outputs
%-----------------------------------------------
%   x_out - array of all x values
%   t_out - array of all t values
%   U_out - array of all solves values over the space of x and time of t
%-----------------------------------------------


%-------Step 0: Argument Checking--------------%   
    
    % kappa should be a scalar 
    if ~isscalar( kappa ) 
        throw( MException( 'MATLAB:invalid_argument', 'the argument kappa is not a scalar' ) );
    end
    
    % x_rng should be a row vector with 2 values    
    if ~all( size( x_rng ) == [1, 2] )
        throw( MException( 'MATLAB:invalid_argument','the argument x_rng is not a row vector with two entries' ) );
    end
    
    % nx should be a scalar    
    if ~isscalar( nx ) 
        throw( MException( 'MATLAB:invalid_argument', 'the argument nx is not a scalar' ) );
    end

    % t_rng should be a row vector with 2 values  
    if ~all( size( t_rng ) == [1, 2] )
        throw( MException( 'MATLAB:invalid_argument','the argument t_rng is not a row vector with two entries' ) );
    end

    % nt should be a scalar     
    if ~isscalar( nt ) 
        throw( MException( 'MATLAB:invalid_argument', 'the argument nt is not a scalar' ) );
    end
    
    % u_init should be a function handle
    if ~isa( u_init, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument','the argument u_init is supposed to be a function handle' ) );
    end

    % u_bndry should be a row vector of function handle
    if ~isa( u_bndry, 'function_handle') 
        throw( MException( 'MATLAB:invalid_argument','the argument u_bndry is supposed to row vector of function handles' ) );
    end

%--------------Step 1: Error Checking-------------%
% To ensure there are no decaying ossilation, kappa*delta T/(h^2) should be
% less than than 0.5

    % creating of h and delta t
    h = (x_rng(2) - x_rng(1))/(nx-1);
    delta_t = (t_rng(2) - t_rng(1))/(nt-1);
    
    % creating x_out and t_out
    x_out = linspace(x_rng(1), x_rng(2), nx);
    t_out = linspace(t_rng(1), t_rng(2), nt);
    
    % creating U (an nx X nt matrix)
    U_out = zeros(nx, nt);
    
    % Calculating r Value
    r = kappa * delta_t / h^2;
    
    % Checking to see if r is resonable
    if r >= 0.5 
        % Calculating the minimum nt needed for r to be resonable 
        nt_minimum = ceil( 2 + kappa * (t_rng(2) - t_rng(1)) / (0.5 * h^2) );
        
        % Throwing a warning 
        warning('MATLAB:questionable_argument','arguments of nt of %f is suboptimum, Use nt of %d', r, nt_minimum);
    end    
%---------------------------------------------%

%-----------Step 2: Initializing--------------------%
    
% Calculating initial conditions 
    u_init_val = u_init(x_out)';
    
    % Calculating the boundary conditions    
    u_bndry_val = u_bndry(t_out(1:end));
    
    % Creating a boundary Matrix
    bndry_matrix = zeros(nx, nt);
    bndry_matrix(1,:) = u_bndry_val(1);
    bndry_matrix(end,:) = u_bndry_val(2);
    
    % Adding all values to U_out matrix
    U_out(:,1) = u_init_val;
    U_out(:,2:end) = U_out(:,2:end) + bndry_matrix(:,2:end);

%---------------------------------------------%  

%-----------Step 3: Solving - --------------------%
    % setting up M matrix
    M_super_diag = diag( -r * ones(nx-3,1), 1 ); 
    M_diag = diag( 2*(r+1) * ones(nx-2,1));
    M_sub_diag = diag( -r * ones(nx-3, 1), -1 );

    % combining all M matrix
    M = M_super_diag + M_diag + M_sub_diag;
    
    % changing M for b boundary condition to be Neumann Condition
    M(end, end-1) = r * -2/3; 
    M(end, end) = 2 + r * 2/3; 

    % resizing boundary value matrix to make dimensions match
    bndry_matrix(2:3,:) = [];
    bndry_matrix(end,:) = 0;
    
    for k = 1:nt-1
        % setting up diff of U as values were added
        U_diff = diff(U_out, 2);
        
        % setting up the b vector 
        b = 2 * U_out(2:end-1, k) + r * U_diff(:,k) + r * bndry_matrix(:,k+1);
        
        % solving for all values at k+1
        U_out(2:end-1, k+1) = M \ b;
               
        % solving for boundary conditoin b
        U_out(end,k+1) = 4/3 * U_out(end-1,k+1) - 1/3  * U_out(end-2, k+1);
    end
%---------------------------------------------------% 
    
%-------------END OF PROGRAM------------------------%
end
