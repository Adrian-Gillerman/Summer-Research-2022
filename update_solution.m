function u = update_solution( Ord,u,un,unm1,sg,cOption,t,dt,ex )
  x = sg.x;
  dx = sg.dx;
  if Ord == 2
    for j = sg.jrange
      up1 = un(j+1);
      u0  = un(j);
      um1 = un(j-1);
      c   = cc(x(j), cOption);
      uxx = (up1 -2*u0 +um1)/dx^2;
      utt = c^2*uxx + ex.f( x(j),t-dt );
      
      u(j) = 2*u0 - unm1(j) + dt^2*utt;
    end
  elseif Ord == 4
    for j = sg.jrange
      up2 = un(j+2);
      up1 = un(j+1);
      u0  = un(j);
      um1 = un(j-1);    
      um2 = un(j-2);
      
      c   = cc(x(j), cOption);
      cp  = ccprime(x(j), cOption);
      cp2 = cc2prime(x(j), cOption);
      
      uxx   = (up1 - 2*u0 + um1)/dx^2;
      uxxx  = (1*up2 - 2*up1 + 2*um1 - 1*um2)/(2*dx^3);
      uxxxx = (1*up2 - 4*up1 + 6*u0  - 4*um1 +1*um2)/dx^4;
      utt   = c^2*(uxx - dx^2/12*uxxxx)...
            + ex.f( x(j),t-dt );
      utttt = c^2*(2*(cp^2 + c*cp2)*uxx...
             + 4*c*cp*uxxx...
             + c^2*uxxxx)...
             + c^2*ex.fxx( x(j),t-dt )...
             + ex.ftt( x(j),t-dt );
      
      u(j) = 2.*u0 - unm1(j)...
           + dt^2*utt ...
           + dt^4/12*utttt;
    end
  else
      error('unsupported order : in spatial_loop')
  end
end
