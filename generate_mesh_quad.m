function [nodes, elements, boundary] = generate_mesh_quad(L, H, Nx, Ny)
% 规则矩形区域四节点网格
% 节点编号：从左下角开始，先沿x方向，再沿y方向
hx = L / Nx; hy = H / Ny;

% 生成节点坐标（确保顺序正确）
nodes = zeros((Nx+1)*(Ny+1), 2);
nid = 0;
for j = 0:Ny  % y方向，从下到上
    for i = 0:Nx  % x方向，从左到右
        nid = nid + 1;
        nodes(nid, 1) = i * hx;
        nodes(nid, 2) = j * hy;
    end
end

% 单元节点号（Q4，逆时针：左下->右下->右上->左上）
% 形函数顺序：节点1(-1,-1), 节点2(1,-1), 节点3(1,1), 节点4(-1,1)
% 节点编号：nid = j*(Nx+1) + i + 1，其中 j=0..Ny, i=0..Nx
elements = zeros(Nx*Ny, 4);
eid = 0;
for j = 0:Ny-1  % 从下到上，j对应单元的下边所在行
    for i = 0:Nx-1  % 从左到右，i对应单元的左边所在列
        % 计算四个角节点的编号
        n1 = j*(Nx+1) + i + 1;           % 左下 (j行, i列)
        n2 = j*(Nx+1) + (i+1) + 1;       % 右下 (j行, i+1列)
        n3 = (j+1)*(Nx+1) + (i+1) + 1;   % 右上 (j+1行, i+1列)
        n4 = (j+1)*(Nx+1) + i + 1;       % 左上 (j+1行, i列)
        eid = eid + 1;
        elements(eid,:) = [n1, n2, n3, n4];
    end
end

% 边界节点集合
tol = 1e-12;
boundary.left   = find(abs(nodes(:,1)) < tol);
boundary.right  = find(abs(nodes(:,1)-L) < tol);
boundary.bottom = find(abs(nodes(:,2)) < tol);
boundary.top    = find(abs(nodes(:,2)-H) < tol);
end

