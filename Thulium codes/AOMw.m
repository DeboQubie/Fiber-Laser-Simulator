function loss = AOMw(t,Li,Lf,tsw)

  if(length(Li)==1)
      Li = Li*ones(length(Lf),1);
  elseif(length(Lf)==1)
      Lf = Lf*ones(length(Li),1);
  end

  w = pi/tsw;                  % omega for sinusoidal variation
  li = 10.^(-Li/10);lf = 10.^(-Lf/10);
  a = (li + lf)/2;
  b = (li - lf)/2;
%   loss = 10*log10((li+lf+cos(w*t)*(li-lf))/2);       % sinusoidal variation
  loss = 10*log10(repmat(a,1,length(t)) + b*cos(w*t));
  loss(:,abs(t)>tsw) = -Lf*ones(1,length(t(abs(t)>tsw)));
end
