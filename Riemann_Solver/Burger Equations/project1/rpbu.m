function ret = rpbu( uL, uR )
usonic = 0.0;        % sonic point
s = 0.5 * (uL + uR); % shock speed from Rankine-Hugoniot
issonic = (uL < 0) & (0 < uR);
ret = issonic * usonic + (~issonic) .* (uL .* (s >= 0) + uR .* (s < 0));
