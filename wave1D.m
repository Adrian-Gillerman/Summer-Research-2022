function [x,u,bc] = wave1D( Ord,Nx,lambda0,tf,cOption,iOption,plottype )
  % Ord:      Order of the scheme
  % Nx:       number of grid points
  % lambda0:  suggested CFL number
  % tf:       final time
  % cOption:  wave speed, 1 for constant, 2 for variable
  % iOption:  1 = sinusoid, 2 = TZ quadratic in space and time
  % plottype: an array where the first element should be
  %            "IC" for just the initial condition
  %            "startfinish" for the initial conditon and final approx
  %            "movie" for every time step
  %            "nothing" for no plotting
  %            and the second element is either
  %            "solution" for the approximation or
  %            "error" for the error to the approximation

  sg     = Spatial_grid( Ord,Nx );
  cmax   = max(cc(sg.x, cOption));
  dt     = lambda0*sg.dx/cmax;
  Nt     = ceil( tf/dt );
  dt     = tf/Nt;
  % lambda = c*dt/sg.dx;
  t      = 0;
  
  %% Set exact solution object
  ex = Exact(cOption, iOption);
  
  %% Set initial conditions over domain interior
  [u,un,unm1] = setICs( Ord,sg,t,dt,ex );

  %% Set BCs
  bc    = BCs( Ord,sg,ex );
  unm1  = bc.setBCs(unm1, t);
  un    = bc.setBCs(un,   t+dt);
  t     = t + dt;

  %% Plot solution at time t = timestep
  uexact = ex.uex( sg.x,t );
  iO = InputOutput( sg,bc,plottype,ex );
  iO.plotIC( un(sg.jrange),uexact(sg.jrange) )
 
  %% Time stepping loop
  for n = 2:Nt
    %% Set time t
    t = t + dt;
      
    %% Loop over interior
    u = update_solution( Ord,u,un,unm1,sg,cOption,t,dt,ex );
       
    %% Set BCs
    u = bc.setBCs(u, t);

    %% Plot to show movie
    uexact  = ex.uex( sg.x,t );
    iO.movie( n,u(sg.jrange),uexact(sg.jrange) );
    
    %% Update old solutions
    unm1     = un;
    un       = u;
  end
  
  %% Plot solution or error at final time 
  iO.plotfinal( u(sg.jrange),uexact(sg.jrange) );
  
  x   = sg.ig;
  u = u(sg.jrange);
  
  return
end
    