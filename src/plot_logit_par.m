%ifig=1;
 fs=15;
 burn_in=(1==1);  it0=1;
 plot_path=(1==1);
  plot_ac=(1==1);
  plot_dichte=(1==1);
 c=['b';'r';'g';'k'];
 c=['r';'b';'g';'k'];
  c=['g';'k';'g';'k'];
 % c=['k';'g';'g';'k'];

  hold_on=(1==1);
  
  for k=2:size(gammc,3)
  thetamc=squeeze(gammc(:,:,k));
  chtheta='\gamma';
 plottheta;
end