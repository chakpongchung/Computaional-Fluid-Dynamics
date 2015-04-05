B=[log(hh');ones(1,length(hh))];
order=(B*B')\(B*log(linearNorm));
% order=order(1);
diff(linearNorm)./diff(hh)
plot(log(hh),log(linearNorm),'o',[-5:1],order(1)*[-5:1]+order(2))