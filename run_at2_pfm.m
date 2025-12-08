% % ** AT2 相场断裂模型（Phase Field Fracture Model）** 
% % ** 基于 MatPFF-main 结构重新构建 **
% % ** 用于模拟裂纹扩展 **
% % ---------------------------------------
% % Last update: 2024
% % Create date: 2024

clear; close all; clc;
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

%% *** 材料参数设置 ***
Para.PFModel = 1; % 1-AT2; 2-AT1; 3-PFCZM
Para.ndim = 2;  % 维度
Para.isStress = 2;  % 1 - 平面应力, 2 - 平面应变
Para.E = 210e3; % 杨氏模量 (N/mm^2)
Para.nu = 0.3;  % 泊松比
Para.lambda = Para.E*Para.nu/((1+Para.nu)*(1-2*Para.nu)); % Lame常数
Para.mu = Para.E/(2*(1+Para.nu)); % Lame常数
Para.Gc = 1.5; % 临界能量释放率 (N/mm) - 减小以加快裂缝扩展
Para.Len = 0.03; % 相场长度尺度 - 增大以让裂缝更明显

%% *** 网格生成 ***
L = 1.0; H = 1.0;
Nx = 40; Ny = 40;
fprintf('生成网格: %d x %d\n', Nx, Ny);
[node, elem, boundary] = generate_mesh_quad(L, H, Nx, Ny);
Para.NNd = size(node, 1);

% 确保节点坐标格式正确（只保留x,y坐标）
if size(node, 2) > 2
    node = node(:, 1:Para.ndim);
end

%% *** 预计算形函数和导数 ***
fprintf('预计算形函数和导数...\n');
[GaussInfo] = shapeFunc_valueDeriv(elem, node, Para);

%% *** 初始相场设置 ***
Phi = zeros(Para.NNd, 1);

% 初始裂缝：左边延伸到边界的垂直裂缝（位于中线附近）
notch_x_max = 0.1*L;  % 裂缝右边界 - 增大初始裂缝长度以更容易扩展
notch_y_center = H/2;  % 裂缝中心y坐标
notch_y_half = 0.003*H;  % 裂缝半高度 - 稍微增大以更明显
notch_ids = find(node(:,1) <= notch_x_max & ...
                 abs(node(:,2) - notch_y_center) <= notch_y_half);
Phi(notch_ids) = 1.0;  % 初始裂缝区域 phi=1

fprintf('初始裂缝：%d 个节点，位置 x <= %.3f, |y-%.3f| <= %.3f\n', ...
    numel(notch_ids), notch_x_max, notch_y_center, notch_y_half);

%% *** 边界条件设置 ***
% 找到左下角节点（用于完全固定）
tol = 1e-12;
bottom_left = find(abs(node(:,1)) < tol & abs(node(:,2)) < tol);
if isempty(bottom_left)
    [~, bottom_left] = min(sum(node.^2, 2));
end

% 设置边界条件：固定底边和左边，允许右边自由变形
fixNode = [];
% 左下角：完全固定
fixNode = [fixNode; bottom_left, 1, 0];  % ux=0
fixNode = [fixNode; bottom_left, 2, 0];  % uy=0
% 左边：ux=0（排除左下角）
left_nodes = setdiff(boundary.left, bottom_left);
fixNode = [fixNode; [left_nodes, ones(numel(left_nodes),1), zeros(numel(left_nodes),1)]];
% 底边：uy=0（排除左下角）
bottom_nodes = setdiff(boundary.bottom, bottom_left);
fixNode = [fixNode; [bottom_nodes, 2*ones(numel(bottom_nodes),1), zeros(numel(bottom_nodes),1)]];
% 右边：不固定，允许自由变形（移除ux=0约束以允许水平变形）
% 顶边：将在载荷步中设置 uy=disp（值>0表示该自由度将被设置为loaddisp）
top_nodes = boundary.top;
fixNode = [fixNode; [top_nodes, 2*ones(numel(top_nodes),1), ones(numel(top_nodes),1)]]; % 值=1表示将被设置为loaddisp

fprintf('边界条件：固定底边和左边，右边自由\n');
fprintf('  左下角节点 %d: ux=0, uy=0（完全固定）\n', bottom_left);
fprintf('  左边: ux=0, 底边: uy=0, 右边: 自由, 顶边: uy=disp（加载）\n');

%% *** 载荷步参数 ***
loadrate = 0.5; % 加载速率 (mm/s) - 增大以加快加载
dt = 1d-4/loadrate; % 初始时间步长
utop = 0.03; % 最大位移 - 显著增大以产生足够的拉伸应变驱动裂缝扩展
nstep = 150; % 载荷步数 - 增加步数以更平滑地加载并观察扩展过程

%% *** 输出目录 ***
outdir = fullfile('AT2_PFM_Code','vtk_output');
if ~exist(outdir,'dir'); mkdir(outdir); end

%% *** 交错迭代求解 ***
fprintf('\n开始交错迭代求解...\n');
AMtol = 1d-4;
loaddisp = 0;

for inc = 1:nstep
    % 调整时间步长 - 在裂缝扩展阶段使用更小的时间步长
    if loaddisp > 0.3*utop
        dt = 5d-6/loadrate;  % 减小时间步长以更精确捕捉裂缝扩展
    end
    
    loaddisp = loaddisp + dt * loadrate;
    if loaddisp > utop
        loaddisp = utop;
    end
    
    % 设置边界条件（动态更新）
    BC = setupBoundaryCondition(fixNode, Para.NNd*Para.ndim, loaddisp);
    
    % 交错迭代
    AMres = 1; it = 0;
    while AMres > AMtol && it < 25
        % 求解位移场
        [Disp] = assembleElasKK(GaussInfo, elem, Phi, Para, BC);
        
        % 更新历史变量
        [GaussInfo, InF] = updateRefEnerg(GaussInfo, elem, Disp, Para);
        
        % 求解相场
        Phiold = Phi;
        Phi = NewtonItPhaseField(GaussInfo, elem, Disp, Phi, Para);
        
        % 确保初始裂缝区域的phi不会减小
        Phi(notch_ids) = max(Phiold(notch_ids), Phi(notch_ids));
        
        % 收敛判据
        AMres = norm(Phi - Phiold) / (norm(Phi) + 1e-10);
        it = it + 1;
    end
    
    % 计算反力
    BDF = sum(InF(BC.BDforce));
    
    % 输出信息
    fprintf('Step %3d: disp=%.4e, Load=%.4e, iter=%d, res=%.2e\n', ...
        inc, loaddisp, full(BDF), it, AMres);
    
    % 输出VTK - 更频繁地输出以观察裂缝扩展过程
    if mod(inc-1, 3) == 0 || inc == nstep
        vtkname = fullfile(outdir, sprintf('step_%03d.vtk', inc));
        write_vtk_at2(vtkname, node, elem, Disp, Phi, []);
        fprintf('  输出VTK: %s\n', vtkname);
    end
    
    if loaddisp >= utop
        break;
    end
end

fprintf('\n计算完成！结果保存在 %s\n', outdir);

