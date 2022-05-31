function speed = cc(x, cOption)
  if cOption == 1
    speed = 0.9;
  elseif cOption == 2
    speed = 1+ 0.5*x.^2;
  else
    error('unsupported cOption : in cc')
  end
end