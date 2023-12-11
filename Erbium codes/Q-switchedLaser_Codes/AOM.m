function loss = AOM(t,Li,Lf,tsw)


  w = pi/tsw;                  % omega for sinusoidal variation
  li = 10^(Li/10);lf = 10^(Lf/10);
  loss = 10*log10((li+lf+(li-lf)*cos(w*t))/2);       % sinusoidal variation
  loss(t>tsw) = Lf;
end
