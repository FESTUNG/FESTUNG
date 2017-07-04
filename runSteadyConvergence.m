function runSteadyConvergence( fac )

problemData.p = 0;
problemData.h = 1;
problemData.isConvergenceRun = true;

pMin = 4;
pMax = 4;
jMin = 0;
jMax = 5;

vecHmax = zeros( jMax+1, 0);
for j=0:jMax
  vecHmax(j+1) = 1 / ( fac * 2^j);
end

err = 1e10 * ones( pMax+1, jMax+1 );

for p=pMin:pMax
  problemData.p = p;
  
  filename = strcat('p', num2str(p), '.txt');
  fileID = fopen(filename,'w');
  writeReportHeader(fileID);
  
  for j=jMin:jMax
    fprintf('Polynomial degree: %d, mesh level: %d \n', p, j );
    problemData.hmax = vecHmax(j+1) ;
    
    resultProblemData = main('hdg_advection_steady', problemData);
    err(p+1, j+1) = resultProblemData.L2error;
    
    fprintf('Polynomial degree: %d\n', p);
    generateCurrentReport( 1, p, j, vecHmax, err );
    generateCurrentReport( fileID, p, j, vecHmax, err );
  end
  fclose(fileID);
end

end


function writeReportHeader( fileID )
%   Write header
fprintf(fileID, 'lvl \tNel \terr \torder\n');
end

function writeReportLine( fileID, j, ne, err, ord )
%   Write line
  fprintf(fileID, '%d \t%d \t%1.10e \t%1.10e \n', j, ne, err, ord );
end

function generateCurrentReport( fileID, p, j, vecHmax, err )
%   Write header
% writeReportHeader(fileID, 'lvl \tNel \terr \torder\n');
ord = 0.;
if j>0
  ord = log(err(p+1, j)/err(p+1, j+1))/log(2) ;
end 
writeReportLine( fileID, j, 1/vecHmax(j+1)*1/vecHmax(j+1)*2, err(p+1, j+1), ord );
end

function generateReport( fileID, p, jMax, vecHmax, err )
%   Write header
writeReportHeader(fileID, 'lvl \tNel \terr \torder\n');
for j=0:jMax
  ord = 0.;
  if j>0
    ord = log(err(p+1, j)/err(p+1, j+1))/log(2) ;
  end 
  writeReportLine( fileID, j, 1/vecHmax(j+1)*1/vecHmax(j+1)*2, err(p+1, j+1), ord );
end 

end
