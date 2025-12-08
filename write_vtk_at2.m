function write_vtk_at2(fname, nodes, elements, u, phi, Hhist)
% 输出 VTK，标量 phi 与向量位移
% 可选：Hhist - 历史变量（如果提供则输出）
if nargin < 6
    Hhist = [];
end

fid = fopen(fname, 'w');
if fid < 0; error('无法写入 %s', fname); end

nn = size(nodes,1);
ne = size(elements,1);
ndof = 2;

fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'AT2 phase-field result\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

fprintf(fid, 'POINTS %d float\n', nn);
fprintf(fid, '%f %f %f\n', [nodes, zeros(nn,1)]');

% 单元
cell_size = 5; % 4 节点 + 1 计数
fprintf(fid, 'CELLS %d %d\n', ne, ne*cell_size);
for e = 1:ne
    conn = elements(e,:) - 1; % VTK 从 0 开始
    fprintf(fid, '4 %d %d %d %d\n', conn);
end

fprintf(fid, 'CELL_TYPES %d\n', ne);
fprintf(fid, '%d\n', 9*ones(ne,1)); % VTK_QUAD

fprintf(fid, 'POINT_DATA %d\n', nn);
fprintf(fid, 'VECTORS displacement float\n');
ux = u(1:ndof:end);
uy = u(2:ndof:end);
fprintf(fid, '%f %f %f\n', [ux, uy, zeros(nn,1)]');

fprintf(fid, 'SCALARS phase float 1\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%f\n', phi);

% 如果提供了历史变量，也输出
if ~isempty(Hhist)
    fprintf(fid, 'SCALARS history float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%f\n', Hhist);
end

fclose(fid);
end

