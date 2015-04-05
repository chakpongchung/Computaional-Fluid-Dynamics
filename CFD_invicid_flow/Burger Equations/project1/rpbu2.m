function ret = rpbu2( uL, uR )

s = 0.5 * (uL + uR);
if uL <= uR,
  if uR <= 0,
    ret = uR;
  else
    
    if uL >= 0,
      ret = uL;
    else
      ret = 0;
    end
  end
else
  if s > 0,
    ret = uL;
  else
    ret = uR;
  end
end