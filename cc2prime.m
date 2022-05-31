function speed = cc2prime(~, cOption)
  if cOption == 1
    speed = 0;
  elseif cOption == 2
    speed = 1;
  else
    error('unsupported cOption : in cc2prime')
  end
end