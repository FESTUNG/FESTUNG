function ret = extrapolateValue(g, triangles, coordinates, values)
N = length(values);
p = (sqrt(8*N+1)-3)/2;
switch p
  case 0
		ret = values;
  case 1
		FinvX1 = 0.5 / g.areaT(triangles) .* ( g.B(triangles,2,2) * coordinates(1) - g.B(triangles,1,2) * coordinates(2) - g.B(triangles,2,2) * g.coordV0T(triangles,1,1) + g.B(triangles,1,2) * g.coordV0T(triangles,1,2) );
		FinvX2 = 0.5 / g.areaT(triangles) .* (-g.B(triangles,2,1) * coordinates(1) + g.B(triangles,1,1) * coordinates(2) + g.B(triangles,2,1) * g.coordV0T(triangles,1,1) - g.B(triangles,1,1) * g.coordV0T(triangles,1,2) );
		ret = values(1) * (1 - FinvX1 - FinvX2) ...
				+ values(2) *      FinvX1           ...
				+ values(3) *               FinvX2;
	otherwise
		error('Local degree of freedom is not supported for extrapolation.')
end % switch
end % function
