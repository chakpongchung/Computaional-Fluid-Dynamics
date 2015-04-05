function ret = uinit( x, ictype )

xshift = 1.0;
switch ictype
 case 1    % shock (uL > uR)
  uL = 1.0;
  uR = 0.5;
  ret = uR + (uL-uR) * ((x-xshift) <= 0.0);
 case 2    % expansion (uL < uR)
  uL = 0.5;
  uR = 1.0;
  ret = uR + (uL-uR) * ((x-xshift) <= 0.0);
 case 3    % sonic expansion (uL < 0 < uR)
  uL =-0.25; 
  uR = 0.5;
  xshift = 4.0;
  ret = uR + (uL-uR) * ((x-xshift) <= 0.0);
 case 4
  uL = 1.0;
  uR = 0.5;
  ret = uR + (uL-uR) * (abs(x-3) <= 1);   %(square pulse)
 case 5

  ret = ( sin(0.25*pi*(x-0)) );  %(sine wave)
 otherwise % cosine (forms shock later)
  ret = 0.5*(1+ cos(pi*(x-xshift))) .* ((x-xshift) <= 1);
end