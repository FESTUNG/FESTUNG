function createHotStart(fileName, dataDisc, tLvl)
%% Open file.
fileName = [fileName, num2str(tLvl), '.txt'];
file     = fopen(fileName, 'wt'); % if this file exists, then overwrite
fprintf(file, '          %.12e\n', dataDisc.');
%% Close file.
fclose(file);
disp(['Hot-Start file ' fileName ' generated.'])
end % function
