function convergencestudy( Ord,h,errarray )
  figure(2)
  lw = 3;
  fs = 16;
  ms = 16;
  
  if     Ord == 2
    color = "g";
  elseif Ord == 4
    color = "b";
  else
    error('unsupported order : in convergencestudy')
  end
  
  loglog( h,errarray,    color+"x",'lineWidth',lw,'MarkerSize',ms )
  hold on
  loglog( h,1.9e0*h.^Ord,color+"-",'lineWidth',lw,'MarkerSize',ms )
  xlabel( 'h' )
  ylabel( 'max error' )
  set(gca,'FontSize',fs)
end


