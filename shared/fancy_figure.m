function fancy_figure()
  
  h = gcf;
  ax = gca;
  
  % fontsize
  %
  propert = findall(h,'-property','FontSize');
  set(propert,'FontSize',20)
  %
  % LATEX stuff
  %
  propert = findall(h,'-property','Interpreter');
  set(propert,'Interpreter','Latex')
  %
  % axis stuff
  %
  ax.TickLength = [0 0];

  % % save stuff
  % %
  % % print( variable, 'name', '-dfileextension', '-rresolution' )
  % %
  % % example
  % %
  % % print(gcf,'figure-example','-dpng','-r650')
  % %
  % print(gcf,'figure-example','-dpdf')

end