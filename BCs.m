classdef BCs
  % BCs Makes BC object that can set the BCs
  %  Ord:      the order of the boundary conditions
  %  sg:       the spatial grid
  %  cOption:  the speed of the wave, either constant or variable
  %  ex:       exact solution object that determines forcing
  %  BCtypes:  the first column determines the "left" or "right" boundary and 
  %            the second column is the dimension.
  %            Either a  "D" for Dirichlet condition or 
  %                   an "N" for Neumann condition
  %
  %  BCs:    constructs a BCs object
  %  setBCs: makes the BCs of u and returns u
    
  properties
    Ord
    sg
    cOption
    ex
    BCtypes = ["D"; "D"];
    BCstr
  end
    
  methods
    function obj = BCs( Ord,sg,ex )
      %BCs Construct an instance of this class
      %  Store the arguments as properties of the class
      obj.Ord     = Ord;
      obj.sg      = sg;
      obj.cOption = ex.cOption; 
      obj.ex      = ex;
      obj.BCstr   = strrep(strjoin(obj.BCtypes),' ','');
    end
        
    function u = setBCs(obj, u, t)
      %setBCs We set the boundary conditions of u at time t
      if obj.BCtypes(1,1) == "D"
        u = setleftDirichlet(obj, u, t);
      elseif obj.BCtypes(1,1) == "N"
        u = setleftNeumann(obj, u, t);
      else
        error('unsupported BC : in BCs')
      end
      if obj.BCtypes(2,1) == "D"
        u = setrightDirichlet(obj, u, t);
      elseif obj.BCtypes(2,1) == "N"
        u = setrightNeumann(obj, u, t);
      else
        error('unsupported BC : in BCs')
      end
    end
    
    function u = setleftDirichlet(obj, u, t)
      dx = obj.sg.dx;
      ja = obj.sg.ja;
      xa = obj.sg.xa;
      exact = obj.ex;
      c   = cc(xa, obj.cOption);
      cp  = ccprime(xa, obj.cOption);
      cp2 = cc2prime(xa, obj.cOption);
      
      if obj.Ord == 2
        u(ja)   = obj.alpha( t );
        u(ja-1) = dx^2/c^2*(obj.alphatt( t )... 
                -                exact.f( xa,t ))...
                + 2*u(ja)...
                - 1*u(ja+1);
      elseif obj.Ord == 4
        A = zeros(5);
        uvec = zeros(5,1);
        uxx   = [0,  1, -2,  1, 0]/dx^2;
        uxxx  = [-1, 2,  0, -2, 1]/(2*dx^3);
        uxxxx = [1, -4,  6, -4, 1]/dx^4;
        a = 2*(cp^2 + c*cp2)*(obj.alphatt(t) - exact.f( xa, t)) +  4*c^3*cp*uxxx;
        % In my code I have a discretization for uxx but in fact we know it
        % exactly
        b = obj.alphatttt(t) - c^2*exact.fxx(xa,t) - exact.ftt(xa,t);
        A(1,:) = c^2*uxx ...
               - dx^2/(12*c^2) ...
               * (b - a);
        A(2,:) = a + c^4*uxxxx;
        A(3,3) = 1;
        A(4,4) = 1;
        A(5,5) = 1;
        uvec(1) = obj.alphatt(t)...
                - exact.f(xa, t);
        uvec(2) = b;
        uvec(3) = obj.alpha(t);
        uvec(4) = u(ja+1);
        uvec(5) = u(ja+2);
        
        unew = A\uvec;
        
        u(ja-2) = unew(1);
        u(ja-1) = unew(2);
        u(ja)   = obj.alpha(t);
      else
        error('unsupported order : in BCs')
      end
    end
    
    function u = setrightDirichlet(obj, u, t)
      dx = obj.sg.dx;
      jb = obj.sg.jb;
      xb = obj.sg.xb;
      exact = obj.ex;
      c   = cc(xb, obj.cOption);
      cp  = ccprime(xb, obj.cOption);
      cp2 = cc2prime(xb, obj.cOption);
      if obj.Ord == 2
        u(jb)   = obj.beta( t );
        u(jb+1) = dx^2/c^2*(obj.betatt( t )... 
                -               exact.f( xb,t ))...
                + 2*u(jb)...
                - 1*u(jb-1);
      elseif obj.Ord == 4
        A    = zeros(5);
        uvec = zeros(5,1);
        B    = zeros(5,1);
        uxx   = [0,  1, -2,  1, 0]/dx^2;
        uxxx  = [-1, 2,  0, -2, 1]/(2*dx^3);
        uxxxx = [1, -4,  6, -4, 1]/dx^4;
        % a = 2*(cp^2 + c*cp2)*(obj.betatt(t) - exact.f( xb, t)) +  4*c^3*cp*uxxx;
        aconst = 2*(cp^2 + c*cp2)*(obj.betatt(t) - exact.f( xb, t));
        avec   = 4*c^3*cp*uxxx;
        % In my code I have a discretization for uxx but in fact we know it
        % exactly
        b = obj.betatttt(t) - c^2*exact.fxx(xb,t) - exact.ftt(xb,t);
        A(1,1) = 1;
        A(2,2) = 1;
        A(3,3) = 1;
        A(4,:) = avec + c^4*uxxxx;
        A(5,:) = c^2*uxx ...
               + dx^2/(12*c^2)*avec;
        uvec(1) = u(jb-2);
        uvec(2) = u(jb-1);
        uvec(3) = obj.beta(t);
        uvec(4) = b;
        uvec(5) = obj.betatt(t);
        B(4)    = aconst;
        B(5)    = exact.f(xb, t)...
                - dx^2/(12*c^2)*(b - aconst);
            
        unew = A\(uvec - B);
        u(jb)   = obj.beta(t);
        u(jb+1) = unew(4);
        % u(jb+1) = exact.uex(xb+dx,t);
        u(jb+2) = unew(5);
      else
        error('unsupported order : in BCs')
      end
    end
    
    function u = setleftNeumann(obj, u, t)
      dx = obj.sg.dx;
      ja = obj.sg.ja;
      xa = obj.sg.xa;
      exact = obj.ex;
      c   = cc(xa, obj.cOption);
      cp  = ccprime(xa, obj.cOption);
      if obj.Ord == 2 
        u(ja-1) = u(ja+1) - 2*dx*obj.Alpha( t );  
      elseif obj.Ord == 4
        uxxx = (obj.Alphatt( t ) - exact.fx( xa,t )...
             - 2*(cp/c)*(obj.alphatt( t ) - exact.f(xa, t)))/c^2;
        u(ja-1) = u(ja+1)...
                - 2*dx*obj.Alpha( t )...
                - dx^3/3*uxxx;
        u(ja-2) = u(ja+2) - 2*u(ja+1) + 2*u(ja-1)...
                - 2*dx^3*uxxx;
      else
        error('unsupported order : in BCs')
      end
    end
    
    function u = setrightNeumann(obj, u, t)
      dx = obj.sg.dx;
      jb = obj.sg.jb;
      xb = obj.sg.xb;
      exact = obj.ex;
      c   = cc(xb, obj.cOption);
      cp  = ccprime(xb, obj.cOption);
      if obj.Ord == 2
        u(jb+1) = u(jb-1) + 2*dx*obj.Beta( t );
      elseif obj.Ord == 4
        uxxx = (obj.Betatt( t ) - exact.fx( xb,t )...
             - 2*(cp/c)*(obj.betatt( t ) - exact.f(xb, t)))/c^2;
        u(jb+1) = u(jb-1) ...
                + 2*dx*obj.Beta( t )...
                + dx^3/3*uxxx;
        u(jb+2) = 2*u(jb+1) - 2*u(jb-1) + u(jb-2)...
                + 2*dx^3*uxxx;
      else
        error('unsupported order : in BCs')
      end
    end
    function z = alpha( obj,t )
      z = obj.ex.uex( obj.sg.xa,t );
    end
    function z = alphatt( obj,t )
      z = obj.ex.utt( obj.sg.xa,t );
    end
    function z = alphatttt( obj,t )
      z = obj.ex.utttt( obj.sg.xa,t );
    end
    function z = Alpha( obj,t )
      z = obj.ex.ux( obj.sg.xa,t );
    end
    function z = Alphatt( obj,t )
      z = obj.ex.uxtt( obj.sg.xa,t );
    end
    function z = beta( obj,t )
      z = obj.ex.uex( obj.sg.xb,t );
    end
    function z = betatt( obj,t )
      z = obj.ex.utt( obj.sg.xb,t );
    end
    function z = betatttt( obj,t )
      z = obj.ex.utttt( obj.sg.xb,t );
    end
    function z = Beta( obj,t )
      z = obj.ex.ux( obj.sg.xb,t );
    end
    function z = Betatt( obj,t )
      z = obj.ex.uxtt( obj.sg.xb,t );
    end
  end
end
