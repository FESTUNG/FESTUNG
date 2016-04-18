refinement = 4;
p = 0:4;
error = zeros(refinement+1,5);
h = 0.3*0.5.^[0 : refinement];
for ord = p
	for k = 0 : refinement
		[error(k+1,1),error(k+1,2),error(k+1,3), error(k+1,4), error(k+1,5)] = mainSWEConv(k, ord);
		order = cell(5,1);
		order{1} = log( error(1:end-1,1) ./ error(2:end,1) ) ./ log( h(1:end-1) ./ h(2:end) );
		order{2} = log( error(1:end-1,2) ./ error(2:end,2) ) ./ log( h(1:end-1) ./ h(2:end) );
		order{3} = log( error(1:end-1,3) ./ error(2:end,3) ) ./ log( h(1:end-1) ./ h(2:end) );
		order{4} = log( error(1:end-1,4) ./ error(2:end,4) ) ./ log( h(1:end-1) ./ h(2:end) );
		order{5} = log( error(1:end-1,5) ./ error(2:end,5) ) ./ log( h(1:end-1) ./ h(2:end) );
		%% clearing all unnecessary variables
		clear h; clear k;
		%% saving error and convergence tables
		save(['p = ' num2str(p)]);
		format long;
		disp(error);
		disp([order{1}, order{2}, order{3}, order{4}, order{5}]);
	end % for
end % for
