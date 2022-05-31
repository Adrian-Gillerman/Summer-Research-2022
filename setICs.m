function [u,un,unm1] = setICs( Ord,sg,t,dt,ex )
  %% Setup u, un, unm1
  NXT     = sg.NXT;
  x       = sg.x;
  dx      = sg.dx;
  cOption = ex.cOption;
  u       = zeros(1, NXT);
  un      = zeros(1, NXT);
  unm1    = zeros(1, NXT);
  if     Ord == 2
    for  j = sg.jrange
      u0m1 = ex.uex( x(j-1),t );
      u00  = ex.uex( x(j),  t );
      u0p1 = ex.uex( x(j+1),t );
      
      c   = cc(x(j), cOption);
      
      u0xx = (u0p1 -2*u00 +u0m1)/dx^2;
      u0tt = c^2*u0xx + ex.f( x(j),t );

      v00  = ex.ut( x(j),t );
    
      unm1(j) = u00;
      un(j)   = u00...
              + dt*v00...
              + dt^2/2*u0tt;
    end
  elseif Ord == 4
    for  j = sg.jrange
      u0m2 = ex.uex( x(j-2),t );
      u0m1 = ex.uex( x(j-1),t );
      u00  = ex.uex( x(j)  ,t );
      u0p1 = ex.uex( x(j+1),t );
      u0p2 = ex.uex( x(j+2),t );
      
      c    = cc(x(j), cOption);
      cp   = ccprime(x(j), cOption);
      cp2  = cc2prime(x(j), cOption);
      
      u0xx   = (-1*u0p2 + 16*u0p1 - 30*u00 + 16*u0m1 - 1*u0m2)/(12*dx^2);
      u0xxx  = ( 1*u0p2 -  2*u0p1 +  2*u0m1 - 1*u0m2)/(2*dx^3);
      u0xxxx = ( 1*u0p2 -  4*u0p1 +  6*u00  - 4*u0m1 + 1*u0m2)/dx^4;
      u0tt   = c^2*u0xx +  ex.f( x(j),t );
      
      u0tttt = c^2*(2*(cp^2 + c*cp2)*u0xx...
             + 4*c*cp*u0xxx...
             + c^2*u0xxxx)...
             + c^2*ex.fxx( x(j),t )...
             + ex.ftt( x(j),t );

      v0m2 = ex.ut( x(j-2),t );
      v0m1 = ex.ut( x(j-1),t );
      v00  = ex.ut( x(j)  ,t );
      v0p1 = ex.ut( x(j+1),t );
      v0p2 = ex.ut( x(j+2),t );
      
      v0xx = (-v0p2+16*v0p1-30*v00+16*v0m1-v0m2)/(12*dx^2);
      v0tt = c^2*v0xx + ex.ft( x(j),t );
    
      unm1(j) = u00;
      un(j)   = u00...
              + dt*v00...
              + dt^2/2 *u0tt...
              + dt^3/6 *v0tt...
              + dt^4/24*u0tttt;
    end
  else
      error('unsupported order : in setICs')
  end
end
