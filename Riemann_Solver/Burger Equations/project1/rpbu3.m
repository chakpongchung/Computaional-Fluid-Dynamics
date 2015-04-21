function ret = rpbu3( uL, uR )
if uL > 0,
  if uR > 0,
    ret = uL;
  else
    s = 0.5 * (uL + uR);
    if s > 0,
      ret = uL;
    else
      ret = uR;
    end
  end
else
  if uR > 0,
    ret = 0;
  else
    ret = uR;
  end
end