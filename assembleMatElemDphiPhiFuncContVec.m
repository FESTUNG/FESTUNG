% TODO
%===============================================================================
%> @file assembleMatElemDphiPhiFuncContVec.m
%>
%> @brief TODO
%===============================================================================
%>
%> @brief TODO
%>
%>
%> 
%> @f[
%> [\mathsf{G}^{m}]_{(k-1)N+i,(k-1)N+j} = \int_{T_{k}}^{} u^m \,\varphi_{kj}\, \partial_{x^{m}} \varphi_{ki} \,\text{d}{\boldsymbol{x}}
%> @f]
%> 
%> @f[
%>   \mathsf{G}^{m}_{T_{K}} \coloneqq
%>   \int_{T_{k}}^{}
%>   u^m
%>   \begin{bmatrix}
%>     \varphi_{k1} \partial_{x^{m}}\, \varphi_{k1} & \dots & \varphi_{kN}\, \partial_{x^{m}} \varphi_{k1}  \\
%>     \vdots & \ddots & \vdots  \\
%>   \varphi_{k1} \,\partial_{x^{m}} \varphi_{kN} & \ldots & \varphi_{kN}\, \partial_{x^{m}} \varphi_{kN}  
%>   \end{bmatrix}
%>   \text{d}{\boldsymbol{x}}\,.
%> @f]
%> 
%> 
%> @f{align*}
%> \int_{T_{k}}^{} u^1(t, {\boldsymbol{x}}) \, \varphi_{kj} \, \partial_{x^{1}} \varphi_{ki} \, \text{d}{\boldsymbol{x}}
%> &\approx
%> \sum_{r=1}^{R} \left(\mathsf{B}^{22}_{k} \, [{\boldsymbol{U}}^1]_{k,r} \, [\mathsf{\hat{G}}]_{1,r,i,j} - \mathsf{B}^{21}_{k} \,[{\boldsymbol{U}}^1]_{k,r}\, [\mathsf{\hat{G}}]_{2,i,j,r}\right) \,, \\
%> \int_{T_{k}}^{} u^2(t, {\boldsymbol{x}}) \, \varphi_{kj} \, \partial_{x^{2}} \varphi_{ki} \, \text{d}{\boldsymbol{x}}
%> &\approx
%> \sum_{r=1}^{R} \left( \mathsf{B}^{11}_{k}\, [{\boldsymbol{U}}^2]_{k,r}\, [\mathsf{\hat{G}}]_{1,r,i,j} - \mathsf{B}^{12}_{k}\, [{\boldsymbol{U}}^2]_{k,r}\, [\mathsf{\hat{G}}]_{2,i,j,r} \right) \,,
%> @f}
%> 
%> 
%> @f[
%> 	[\mathsf{\hat{G}}]_{m,i,j,r} \coloneqq \omega_{r} \, \partial_{\hat{x}^{m}}{} \, \hat{\varphi}_{ki} \, \hat{\varphi}_{kj} \,,
%>   \qquad
%>   [{\boldsymbol{U}}^m]_{k,r} \coloneqq u^m(t, {\boldsymbol{F}}_k({\boldsymbol{q}}_r))\,.
%> @f]
%> 
%> 
%> @f{align*}
%> \mathsf{G}^{1}_{{T_{k}}} 
%> &= 
%> \sum_{r=1}^{R} 
%> \left(
%> B^{22}_{k}
%> \begin{bmatrix}
%> 	\omega_{r} \, u^1(t,{\boldsymbol{F}}_k({\boldsymbol{q}}_r)) \, \partial_{x^{1}} \hat{\varphi}_{1}\, \hat{\varphi}_{1}  & \cdots & \omega_{r} \,\, u^1(t,{\boldsymbol{F}}_k({\boldsymbol{q}}_r)) \, \partial_{x^{1}} \hat{\varphi}_{1}\, \hat{\varphi}_{N} 
%> 	\\
%> 	\vdots  & \ddots & \vdots 
%> 	\\
%> 	\omega_{r} \, u^1(t,{\boldsymbol{F}}_k({\boldsymbol{q}}_r)) \,   \partial_{x^{1}} \hat{\varphi}_{N}\, \hat{\varphi}_{1}  & \cdots & \omega_{r} \, u^1(t,{\boldsymbol{F}}_k({\boldsymbol{q}}_r)) \,  \partial_{x^{1}} \hat{\varphi}_{N}\, \hat{\varphi}_{N} 
%> \end{bmatrix}\right.
%> \\
%> &\phantom{= \sum_{r=1}^{R} }
%> \left.
%> - B^{21}_{k}
%> \begin{bmatrix}
%> \omega_{r} \, u^1(t,{\boldsymbol{F}}_k({\boldsymbol{q}}_r)) \, \partial_{x^{2}} \hat{\varphi}_{1} \, \hat{\varphi}_{1}  & \cdots & \omega_{r} \, u^1(t,{\boldsymbol{F}}_k({\boldsymbol{q}}_r)) \,  \partial_{x^{2}} \hat{\varphi}_{1} \,\hat{\varphi}_{N} 
%> \\
%> \vdots  & \ddots & \vdots 
%> \\
%> \omega_{r} \, u^1(t,{\boldsymbol{F}}_k({\boldsymbol{q}}_r)) \, \partial_{x^{2}} \hat{\varphi}_{N}\,\hat{\varphi}_{1}  & \cdots & \omega_{r} \, u^1(t,{\boldsymbol{F}}_k({\boldsymbol{q}}_r)) \,\partial_{x^{2}} \hat{\varphi}_{N} \,\hat{\varphi}_{N} 
%> \end{bmatrix}
%> \right)
%> \\
%> &= 
%> \sum_{r=1}^{R} \left( B^{22}_{k} \,[{\boldsymbol{U}}^1]_{k,r} \,[\mathsf{\hat{G}}]_{1,:,:,r} - B^{21}_{k}\, [{\boldsymbol{U}}^1]_{k,r} \,[\mathsf{\hat{G}}]_{2,:,:,r} \right)\,.
%> @f}
%> 
%> asdf
%> @f[
%> \mathsf{G}^1 = 
%> \sum_{r=1}^{R} \mathsf{I}_{K\times K} \otimes_\mathrm{V} \left( 
%>   \begin{bmatrix} B_1^{22} \, [{\boldsymbol{U}}^1]_{1,r} \\ \vdots \\ B_K^{22} [{\boldsymbol{U}}^1]_{K,r} \end{bmatrix}  \otimes [\mathsf{\hat{G}}]_{1,:,:,r} -
%>   \begin{bmatrix} B_1^{21} \, [{\boldsymbol{U}}^1]_{1,r} \\ \vdots \\ B_K^{21} [{\boldsymbol{U}}^1]_{K,r} \end{bmatrix}  \otimes [\mathsf{\hat{G}}]_{2,:,:,r}
%> \right) \,,
%> @f]
%>
%>
%>
%> 
%> All other entries are zero.
%> @param  TODO TODO
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Alexander Jaust, 2017
%> @author Balthasar Reuter, 2017
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
function ret = assembleMatElemDphiPhiFuncContVec(g, refElemDphiPhiPerQuad, funcCont1, funcCont2, qOrd)
K = g.numT;
[N, ~, R] = size(refElemDphiPhiPerQuad{1});
validateattributes(refElemDphiPhiPerQuad, {'cell'}, {'size', [2 1]}, mfilename, 'refElemDphiPhiPerQuad');
validateattributes(refElemDphiPhiPerQuad{1}, {'numeric'}, {'size', [N N R]}, mfilename, 'refElemDphiPhiPerQuad{1}');
validateattributes(refElemDphiPhiPerQuad{2}, {'numeric'}, {'size', [N N R]}, mfilename, 'refElemDphiPhiPerQuad{2}');
validateattributes(funcCont1, {'function_handle'}, {'size', [1 1]}, mfilename, 'funcCont1');
validateattributes(funcCont2, {'function_handle'}, {'size', [1 1]}, mfilename, 'funcCont2');

if nargin < 5, p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1); end
[Q1, Q2, ~] = quadRule2D(qOrd);

% Assemble matrix
ret = { zeros(K*N, N), zeros(K*N,N) };
for r = 1 : R
  valOnQuad1 = funcCont1(g.mapRef2Phy(1, Q1(r), Q2(r)), g.mapRef2Phy(2, Q1(r), Q2(r)));
  ret{1} = ret{1} + kron(g.B(:,2,2) .* valOnQuad1, refElemDphiPhiPerQuad{1}(:, :, r)) ...
                  - kron(g.B(:,2,1) .* valOnQuad1, refElemDphiPhiPerQuad{2}(:, :, r));
        
  valOnQuad2 = funcCont2(g.mapRef2Phy(1, Q1(r), Q2(r)), g.mapRef2Phy(2, Q1(r), Q2(r)));
  ret{2} = ret{2} - kron(g.B(:,1,2) .* valOnQuad2, refElemDphiPhiPerQuad{1}(:, :, r)) ...
                  + kron(g.B(:,1,1) .* valOnQuad2, refElemDphiPhiPerQuad{2}(:, :, r));
end %for
ret{1} = kronVec(speye(K,K), ret{1});
ret{2} = kronVec(speye(K,K), ret{2});
end
