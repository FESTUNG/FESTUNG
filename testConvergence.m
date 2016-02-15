lvl = 0:6;
p = 0;%1:4;
typeSlopeLim = cell(4,1);
typeSlopeLim{1} = 'linear';
typeSlopeLim{2} = 'kuzmin';
typeSlopeLim{3} = 'strict';
typeSlopeLim{4} = 'edge';

for ord = p
  err = zeros(size(lvl));

  for i = 1:length(lvl)
    err(i) = mainAdvectionConv(2.^(-lvl(i)), ord, false, 'none');
  end

  eoc = log( err(1:end-1) ./ err(2:end) ) ./ log( 2.^(-lvl(1:end-1)) ./ (2.^(-lvl(2:end))) );

  f = fopen(['p' num2str(ord) '.txt'],'w');
  fprintf(f, 'lvl   error      eoc\n');
  fprintf(f, '%3d  %10.3e  %4.2f\n', [lvl; err; [0 eoc]]);
  fclose(f);
end
return
for j = 1:size(typeSlopeLim,1)
  type = typeSlopeLim{j};
  for ord = p
    err = zeros(size(lvl));

    for i = 1:length(lvl)
      err(i) = mainAdvectionConv(2.^(-lvl(i)), ord, true, type);
    end

    eoc = log( err(1:end-1) ./ err(2:end) ) ./ log( 2.^(-lvl(1:end-1)) ./ (2.^(-lvl(2:end))) );

    f = fopen([type '_p' num2str(ord) '.txt'],'w');
    fprintf(f, 'lvl   error      eoc\n');
    fprintf(f, '%3d  %10.3e  %4.2f\n', [lvl; err; [0 eoc]]);
    fclose(f);
  end
end