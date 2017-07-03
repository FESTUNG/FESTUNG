function runSteadyConvergence

problemData.p = 0;
problemData.h = 1;

pMax = 4;
jMax = 7;

vecHmax = zeros( jMax+1, 0);
for j=0:jMax
  vecHmax(j+1) = 1 / ( 3 * 2^j);
end

err = 1e10 * ones( pMax+1, jMax+1 );

for p=0:pMax
  problemData.p = p;
  for j=0:jMax
    fprintf('Polynomial degree: %d, mesh level: %d \n', p, j );
    problemData.hmax = vecHmax(j+1) ;
    
    resultProblemData = main('hdg_advection_steady', problemData);
    err(p+1, j+1) = resultProblemData.L2error;
  end

  fprintf('Polynomial degree: %d\n', p);
  generateReport( 1, p, jMax, vecHmax, err );
  
  filename = strcat('p', num2str(p), '.txt');
  fileID = fopen(filename,'w');
  generateReport( fileID, p, jMax, vecHmax, err );
  fclose(fileID);
end

end

function generateReport( fileID, p, jMax, vecHmax, err )


%   Write header
fprintf(fileID, 'lvl \tNel \terr \torder\n');
for j=0:jMax
  if j>0
    fprintf(fileID, '%d \t%d \t%g \t%g \n', j, 1/vecHmax(j+1)*1/vecHmax(j+1)*2, err(p+1, j+1), log(err(p+1, j)/err(p+1, j+1))/log(2) );
  else
    fprintf(fileID, '%d \t%d \t%g \t0.0 \n', j, 1/vecHmax(j+1)*1/vecHmax(j+1)*2, err(p+1, j+1) );
  end
end 

end
