function BoundaryCondition = setupBoundaryCondition(fixNode, NDof, loaddisp)
% 设置边界条件结构
% 输入: fixNode - [节点号, 自由度(1-x, 2-y), 值]（值>0表示该自由度将被设置为loaddisp）
%       NDof - 总自由度数
%       loaddisp - 当前载荷步的位移值
% 输出: BoundaryCondition - 边界条件结构

% 更新有非零值的边界节点
fixNode(fixNode(:,3)>0, 3) = loaddisp;
MoveNode = fixNode(fixNode(:,3)>0, 1:2);

% Neumann(Natural) Boundary Condition
F = zeros(NDof, 1);

% BCs
BoundaryCondition.DirchletDOF = fixNode(:,1)*2 + (fixNode(:,2)-2);
BoundaryCondition.Dirichlet   = fixNode(:,3);
BoundaryCondition.FreeDOF     = setdiff([1:NDof]', BoundaryCondition.DirchletDOF);
BoundaryCondition.BDforce     = MoveNode(:,1)*2 + (MoveNode(:,2)-2);
BoundaryCondition.RHS         = F;

end

