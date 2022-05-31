function speed = ccprime(x, cOption)
  if cOption == 1
    speed = 0;
  elseif cOption == 2
    speed = x;
  else
    error('unsupported cOption : in ccprime')
  end
end