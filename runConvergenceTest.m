% Script for checking the convergence order for the different unknowns of 
% the shllow water equations.
%
%===============================================================================
%> @file runConvergenceTest.m
%>
%> @brief Script for checking the convergence order for the different
%>        unknowns of the shllow water equations.
%===============================================================================
%>
%> @brief Script for checking the convergence order for the different
%>        unknowns of the shllow water equations.
%>
%> A structured grid is created. For the purpose of bandwidth it makes 
%> sense to have fewer elements in y- than in x-direction. 
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2018
%> 
%> @par License
%> @parblock
%> This program is free software: you can redistribute it and/or modify
%> it under the terms of the GNU General Public License as published by
%> the Free Software Foundation, either version 3 of the License, or
%> (at your option) any later version.
%>
%> This program is distributed in the hope that it will be useful,
%> but WITHOUT ANY WARRANTY; without even the implied warranty of
%> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%> GNU General Public License for more details.
%>
%> You should have received a copy of the GNU General Public License
%> along with this program.  If not, see <http://www.gnu.org/licenses/>.
%> @endparblock
%
clc
clear
cd ..
maxLvl = 4;
hmax = 10;
h = hmax * 0.5.^(1:maxLvl)';
table = zeros(maxLvl+1, 10);
for lvl = 0:maxLvl
  pd.refinement = lvl;
  out = main('swe', pd);
  table(lvl+1, 1:2:end-1) = out.err([1 3:end]);
  if lvl > 0
    table(lvl+1,2:2:end) = log(table(lvl, 1:2:end-1) ./ table(lvl+1, 1:2:end-1)) / log(2);
  end % if
end % for
fileID = fopen(['convergence_results_' datestr(now,'yyyy-mm-dd HH:MM:SS') '.txt'], 'wt');
fprintf(fileID, 'Errors:\n   $h$ &      $e_H$ &    EOC$_H$ &   $e_{uH}$ & EOC$_{uH}$ &   $e_{vH}$ & EOC$_{vH}$ &      $e_u$ &    EOC$_u$ &      $e_v$ &    EOC$_v$ \\\\ \n');
fprintf(fileID, '-----------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf(fileID, '%1.3f &   %1.2E &         -- &   %1.2E &         -- &   %1.2E &         -- &   %1.2E &         -- &   %1.2E &         -- \\\\ \n', [hmax table(1,1:2:end-1)]);
fprintf(fileID, ' %1.3f &   %1.2E &       %1.2f &   %1.2E &       %1.2f &   %1.2E &       %1.2f &   %1.2E &       %1.2f &   %1.2E &       %1.2f \\\\ \n', [h table(2:end,:)]');
fclose(fileID);
cd swe
