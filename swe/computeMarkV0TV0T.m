% Computes the cell markV0TV0T for a structured grid grid over a 
% rectangular domain.
%
%===============================================================================
%> @file
%>
%> @brief Computes the cell markV0TV0T for a structured grid grid over a 
%>        rectangular domain using a slow but memory efficient algorithm.
%===============================================================================
%>
%> @brief Computes the cell markV0TV0T for a structured grid grid over a 
%>        rectangular domain using a slow but memory efficient algorithm.
%>
%> For the special case of a structured triangulation grid over a
%> rectangular domain, the sparse matrices in markV0TV0T can be computed 
%> without a large meomory requirement but in a very slow iterative and 
%> process and by adressing certain entries of sparse matrices. 
%> Instead of calling this routine repeatedly, it is wiser to save grid
%> data due to the otherwise significantly longer preprocessing time.
%>
%> @param  numTx        The number of triangles in x-direction.
%> @param  numTy        The number of triangles in y-direction.
%> @param  isPeriodic   Logical that is used to decide whether to impose
%>                      periodic boundary conditions.
%> @retval markV0TV0T   The @f$(n^-, n^+)@f$th entry of this cell is a 
%>                      sparse @f$K \times K@f$ array whose @f$(k^-,k^+)@f$th
%>                      entry is one if @f$v_{k^-n^-}=v_{k^+n^+}@f$
%>                      @f$[3 \times 3 \text{ (cell)}]@f$.
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
function markV0TV0T = computeMarkV0TV0T(numTx, numTy, isPeriodic)

if nargin < 3
  isPeriodic = false;
end
if nargin < 2
  numTy = numTx;
end

K = 2*numTx*numTy;
markV0TV0T = cell(3,3);
markV0TV0T{1,1} = sparse(K,K);
markV0TV0T{1,2} = sparse(K,K);
markV0TV0T{1,3} = sparse(K,K);
markV0TV0T{2,1} = sparse(K,K);
markV0TV0T{2,2} = sparse(K,K);
markV0TV0T{2,3} = sparse(K,K);
markV0TV0T{3,1} = sparse(K,K);
markV0TV0T{3,2} = sparse(K,K);
markV0TV0T{3,3} = sparse(K,K);

for kx = 1 : numTx-1
  for ky = 2 : 2 : 2*numTy-2
    
    k = (kx-1)*2*numTy + ky;
    
    markV0TV0T{1,1}(k+2,k+2) = 1;
    markV0TV0T{1,1}(k+2,k+2*numTy+1) = 1;
    markV0TV0T{1,1}(k+2*numTy+1,k+2) = 1;
    markV0TV0T{1,1}(k+2*numTy+1,k+2*numTy+1) = 1;
    
    markV0TV0T{1,2}(k+2,k+1) = 1;
    markV0TV0T{1,2}(k+2,k) = 1;
    markV0TV0T{1,2}(k+2*numTy+1,k+1) = 1;
    markV0TV0T{1,2}(k+2*numTy+1,k) = 1;
    
    markV0TV0T{1,3}(k+2,k+2*numTy-1) = 1;
    markV0TV0T{1,3}(k+2,k+2*numTy) = 1;
    markV0TV0T{1,3}(k+2*numTy+1,k+2*numTy-1) = 1;
    markV0TV0T{1,3}(k+2*numTy+1,k+2*numTy) = 1;
    
    markV0TV0T{2,1}(k,k+2) = 1;
    markV0TV0T{2,1}(k,k+2*numTy+1) = 1;
    markV0TV0T{2,1}(k+1,k+2) = 1;
    markV0TV0T{2,1}(k+1,k+2*numTy+1) = 1;
    
    markV0TV0T{2,2}(k,k) = 1;
    markV0TV0T{2,2}(k,k+1) = 1;
    markV0TV0T{2,2}(k+1,k) = 1;
    markV0TV0T{2,2}(k+1,k+1) = 1;
    
    markV0TV0T{2,3}(k,k+2*numTy-1) = 1;
    markV0TV0T{2,3}(k,k+2*numTy) = 1;
    markV0TV0T{2,3}(k+1,k+2*numTy-1) = 1;
    markV0TV0T{2,3}(k+1,k+2*numTy) = 1;
    
    markV0TV0T{3,1}(k+2*numTy-1,k+2) = 1;
    markV0TV0T{3,1}(k+2*numTy-1,k+2*numTy+1) = 1;
    markV0TV0T{3,1}(k+2*numTy,k+2) = 1;
    markV0TV0T{3,1}(k+2*numTy,k+2*numTy+1) = 1;
    
    markV0TV0T{3,2}(k+2*numTy-1,k) = 1;
    markV0TV0T{3,2}(k+2*numTy-1,k+1) = 1;
    markV0TV0T{3,2}(k+2*numTy,k) = 1;
    markV0TV0T{3,2}(k+2*numTy,k+1) = 1;
    
    markV0TV0T{3,3}(k+2*numTy-1,k+2*numTy-1) = 1;
    markV0TV0T{3,3}(k+2*numTy,k+2*numTy) = 1;
    markV0TV0T{3,3}(k+2*numTy-1,k+2*numTy) = 1;
    markV0TV0T{3,3}(k+2*numTy,k+2*numTy-1) = 1;
    
  end
end

for k = 1 : 2*numTy : 2*numTy*(numTx-2)+1 % south bdr
  markV0TV0T{1,1}(k+1,k+1) = 1;
  markV0TV0T{1,1}(k+1,k+2*numTy) = 1;
  markV0TV0T{1,1}(k+2*numTy,k+1) = 1;
  markV0TV0T{1,1}(k+2*numTy,k+2*numTy) = 1;
  
  markV0TV0T{1,2}(k+1,k) = 1;
  markV0TV0T{1,2}(k+2*numTy,k) = 1;
  
  markV0TV0T{2,1}(k,k+1) = 1;
  markV0TV0T{2,1}(k,k+2*numTy) = 1;
  
  markV0TV0T{2,2}(k,k) = 1;
  
  if isPeriodic
    markV0TV0T{1,2}(k+1,k+2*numTy-1) = 1;
    markV0TV0T{1,2}(k+2*numTy,k+2*numTy-1) = 1;
    
    markV0TV0T{1,3}(k+1,k+4*numTy-2) = 1;
    markV0TV0T{1,3}(k+1,k+4*numTy-1) = 1;
    markV0TV0T{1,3}(k+2*numTy,k+4*numTy-2) = 1;
    markV0TV0T{1,3}(k+2*numTy,k+4*numTy-1) = 1;
    
    markV0TV0T{2,1}(k+2*numTy-1,k+1) = 1;
    markV0TV0T{2,1}(k+2*numTy-1,k+2*numTy) = 1;
    
    markV0TV0T{2,2}(k,k+2*numTy-1) = 1;
    markV0TV0T{2,2}(k+2*numTy-1,k) = 1;
    
    markV0TV0T{2,3}(k,k+4*numTy-2) = 1;
    markV0TV0T{2,3}(k,k+4*numTy-1) = 1;
    
    markV0TV0T{3,1}(k+4*numTy-2,k+1) = 1;
    markV0TV0T{3,1}(k+4*numTy-2,k+2*numTy) = 1;
    markV0TV0T{3,1}(k+4*numTy-1,k+1) = 1;
    markV0TV0T{3,1}(k+4*numTy-1,k+2*numTy) = 1;
    
    markV0TV0T{3,2}(k+4*numTy-2,k) = 1;
    markV0TV0T{3,2}(k+4*numTy-1,k) = 1;
  end
end

for k = 2*numTy*(numTx-1)+2 : 2 : 2*numTy*numTx-2 % east bdr
  markV0TV0T{1,1}(k+2,k+2) = 1;
  
  markV0TV0T{1,2}(k+2,k+1) = 1;
  markV0TV0T{1,2}(k+2,k) = 1;
  
  markV0TV0T{2,1}(k,k+2) = 1;
  markV0TV0T{2,1}(k+1,k+2) = 1;
  
  markV0TV0T{2,2}(k,k) = 1;
  markV0TV0T{2,2}(k+1,k) = 1;
  markV0TV0T{2,2}(k,k+1) = 1;
  markV0TV0T{2,2}(k+1,k+1) = 1;
end

for k = 2*numTy : 2*numTy : 2*numTy*(numTx-1) % north bdr
  markV0TV0T{2,2}(k,k) = 1;
  
  markV0TV0T{2,3}(k,k+2*numTy-1) = 1;
  markV0TV0T{2,3}(k,k+2*numTy) = 1;
  
  markV0TV0T{3,2}(k+2*numTy-1,k) = 1;
  markV0TV0T{3,2}(k+2*numTy,k) = 1;
  
  markV0TV0T{3,3}(k+2*numTy-1,k+2*numTy-1) = 1;
  markV0TV0T{3,3}(k+2*numTy,k+2*numTy-1) = 1;
  markV0TV0T{3,3}(k+2*numTy-1,k+2*numTy) = 1;
  markV0TV0T{3,3}(k+2*numTy,k+2*numTy) = 1;
end

for k = 1 : 2 : 2*numTy-3 % west bdr
  markV0TV0T{1,1}(k+2,k+2) = 1;
  
  markV0TV0T{1,3}(k+2,k+1) = 1;
  markV0TV0T{1,3}(k+2,k) = 1;
  
  markV0TV0T{3,1}(k+1,k+2) = 1;
  markV0TV0T{3,1}(k,k+2) = 1;
  
  markV0TV0T{3,3}(k,k) = 1;
  markV0TV0T{3,3}(k,k+1) = 1;
  markV0TV0T{3,3}(k+1,k) = 1;
  markV0TV0T{3,3}(k+1,k+1) = 1;
  
  if isPeriodic
    markV0TV0T{1,1}(k+2,k+2*numTy*(numTx-1)+3) = 1;
    markV0TV0T{1,1}(k+2*numTy*(numTx-1)+3,k+2) = 1;
    
    markV0TV0T{1,2}(k+2,k+2*numTy*(numTx-1)+1) = 1;
    markV0TV0T{1,2}(k+2,k+2*numTy*(numTx-1)+2) = 1;
    
    markV0TV0T{1,3}(k+2*numTy*(numTx-1)+3,k) = 1;
    markV0TV0T{1,3}(k+2*numTy*(numTx-1)+3,k+1) = 1;
    
    markV0TV0T{2,1}(k+2*numTy*(numTx-1)+2,k+2) = 1;
    markV0TV0T{2,1}(k+2*numTy*(numTx-1)+1,k+2) = 1;
    
    markV0TV0T{2,3}(k+2*numTy*(numTx-1)+2,k) = 1;
    markV0TV0T{2,3}(k+2*numTy*(numTx-1)+1,k) = 1;
    markV0TV0T{2,3}(k+2*numTy*(numTx-1)+2,k+1) = 1;
    markV0TV0T{2,3}(k+2*numTy*(numTx-1)+1,k+1) = 1;
    
    markV0TV0T{3,1}(k,k+2*numTy*(numTx-1)+3) = 1;
    markV0TV0T{3,1}(k+1,k+2*numTy*(numTx-1)+3) = 1;
    
    markV0TV0T{3,2}(k,k+2*numTy*(numTx-1)+1) = 1;
    markV0TV0T{3,2}(k,k+2*numTy*(numTx-1)+2) = 1;
    markV0TV0T{3,2}(k+1,k+2*numTy*(numTx-1)+1) = 1;
    markV0TV0T{3,2}(k+1,k+2*numTy*(numTx-1)+2) = 1;
  end
end


% south east corner
markV0TV0T{1,1}(2*numTy*(numTx-1)+2,2*numTy*(numTx-1)+2) = 1;
markV0TV0T{1,2}(2*numTy*(numTx-1)+2,2*numTy*(numTx-1)+1) = 1;
markV0TV0T{2,1}(2*numTy*(numTx-1)+1,2*numTy*(numTx-1)+2) = 1;
markV0TV0T{2,2}(2*numTy*(numTx-1)+1,2*numTy*(numTx-1)+1) = 1;

% north east corner
markV0TV0T{2,2}(2*numTy*numTx,2*numTy*numTx) = 1;

% north west corner
markV0TV0T{3,3}(2*numTy,2*numTy) = 1;
markV0TV0T{3,3}(2*numTy-1,2*numTy) = 1;
markV0TV0T{3,3}(2*numTy,2*numTy-1) = 1;
markV0TV0T{3,3}(2*numTy-1,2*numTy-1) = 1;

% south west corner
markV0TV0T{1,1}(1,1) = 1;

if isPeriodic
  markV0TV0T{1,1}(1,2+2*numTy*(numTx-1)) = 1;
  markV0TV0T{1,1}(2+2*numTy*(numTx-1),1) = 1;
  
  markV0TV0T{1,2}(1,1+2*numTy*(numTx-1)) = 1;
  markV0TV0T{1,2}(1,2*numTy*numTx) = 1;
  markV0TV0T{1,2}(2+2*numTy*(numTx-1),2*numTx*numTy) = 1;
  
  markV0TV0T{1,3}(1,2*numTy-1) = 1;
  markV0TV0T{1,3}(1,2*numTy) = 1;
  markV0TV0T{1,3}(2*numTy*(numTx-1)+2,2*numTy-1) = 1;
  markV0TV0T{1,3}(2*numTy*(numTx-1)+2,2*numTy) = 1;
  
  markV0TV0T{2,1}(2*numTy*(numTx-1)+1,1) = 1;
  markV0TV0T{2,1}(2*numTy*numTx,1) = 1;
  markV0TV0T{2,1}(2*numTy*numTx,2*numTy*(numTx-1)+2) = 1;
  
  markV0TV0T{2,2}(2*numTy*(numTx-1)+1,2*numTy*numTx) = 1;
  markV0TV0T{2,2}(2*numTy*numTx,2*numTy*(numTx-1)+1) = 1;
  
  markV0TV0T{2,3}(2*numTy*(numTx-1)+1,2*numTy-1) = 1;
  markV0TV0T{2,3}(2*numTy*(numTx-1)+1,2*numTy) = 1;
  markV0TV0T{2,3}(2*numTy*numTx,2*numTy-1) = 1;
  markV0TV0T{2,3}(2*numTy*numTx,2*numTy) = 1;
  
  markV0TV0T{3,1}(2*numTy-1,1) = 1;
  markV0TV0T{3,1}(2*numTy,1) = 1;
  markV0TV0T{3,1}(2*numTy-1,2*numTy*(numTx-1)+2) = 1;
  markV0TV0T{3,1}(2*numTy,2*numTy*(numTx-1)+2) = 1;
  
  markV0TV0T{3,2}(2*numTy-1,2*numTy*(numTx-1)+1) = 1;
  markV0TV0T{3,2}(2*numTy,2*numTy*(numTx-1)+1) = 1;
  markV0TV0T{3,2}(2*numTy-1,2*numTy*numTx) = 1;
  markV0TV0T{3,2}(2*numTy,2*numTy*numTx) = 1;
end


markV0TV0T{1,1} = logical(markV0TV0T{1,1});
markV0TV0T{1,2} = logical(markV0TV0T{1,2});
markV0TV0T{1,3} = logical(markV0TV0T{1,3});
markV0TV0T{2,1} = logical(markV0TV0T{2,1});
markV0TV0T{2,2} = logical(markV0TV0T{2,2});
markV0TV0T{2,3} = logical(markV0TV0T{2,3});
markV0TV0T{3,1} = logical(markV0TV0T{3,1});
markV0TV0T{3,2} = logical(markV0TV0T{3,2});
markV0TV0T{3,3} = logical(markV0TV0T{3,3});

if max(numTx,numTy) < 50
  testThis(markV0TV0T, numTx, numTy);
end

end % function

function testThis(ret, nx, ny)

g = domainRectanglePeriodic(0:1,0:1,1/nx,1/ny);
for nn = 1:3
  for np = 1:3
    if max(max(abs(g.markV0TV0T{nn,np} - ret{nn,np}))) ~= 0
      error(['error at (nn,np) = (', num2str(nn) ',' num2str(np) ')'])
    end
  end
end

end
