addpath(genpath("../Graphs"));

% Set some test cases with known outcomes.

% 2D: Nearest neighbour in x, nearest neighbour in y, nearest neighbours,
% diagonal neighbours, isolated, with both PBC and OBC.

G_2d_nx_pb = HypCub([4 4],[1 1],[1 0],1);
G_2d_nx_oxb = HypCub([4 4],[0 1],[1 0],1);
G_2d_ny_pb = HypCub([4 4],[1 1],[0 1],1);
G_2d_ny_oyb = HypCub([4 4],[0 0],[0 1],1);
G_2d_nn_pb = HypCub([4 4],[1 1],eye(2),1);
G_2d_nn_ob = HypCub([4 4],[0 0],eye(2),1);
G_2d_dnx_pb = HypCub([4 4],[1 1],[1, -1; 1, 1],1);
G_2d_dnx_ob = HypCub([4 4],[0 0],[1, -1; 1, 1],1);
G_2d_dny_pb = HypCub([4 4],[1 1],[1, 1; -1, 1],1);
G_2d_dny_ob = HypCub([4 4],[0 0],[1, 1; -1, 1],1);

G_2d_iso = HypCub([4 4],[0 0],[0 0],1);

GraphArray = {G_2d_nx_pb; G_2d_nx_oxb; G_2d_ny_pb; G_2d_ny_oyb; G_2d_nn_pb; ...
    G_2d_nn_ob; G_2d_dnx_pb; G_2d_dnx_ob; G_2d_dny_pb; G_2d_dny_ob; G_2d_iso};

GraphIDs = {'2D x-neighbours, PBC'; '2D x-neighbours, OBC'; '2D y-neighbours, PBC'; ...
    '2D y-neighbours, OBC'; '2D nearest neighbours, PBC'; '2D nearest neighbours, OBC'; ...
    '2D diagonal neighbours, positive-x, PBC'; '2D diagonal neighbours, positive-x, OBC'; ...
    '2D diagonal neighbours, positive-y, PBC'; '2D diagonal neighbours, positive-y, OBC'; ...
    '2D isolated sites'};

BoundarySites = [1, 4, 13, 16];
EmptyBonds = [0, 0, 0, 0;
              0, 1, 0, 1;
              0, 0, 0, 0;
              0, 0, 1, 1;
              0, 0, 0, 0;
              0, 1, 1, 2;
              0, 0, 0, 0;
              1, 2, 1, 2;
              0, 0, 0, 0;
              1, 1, 2, 2;
              0, 0, 0, 0];

SubLatNums = [4; 4; 4; 4; 1; 1; 2; 2; 2; 2; 16];

EmptyBondNums = [0; 4; 0; 4; 0; 8; 0; 14; 0; 14; 0];

%% Test 1: Edge of lattice neighbours

for g = 1:numel(GraphArray)
    disp(['Testing Graph: ' GraphIDs{g}]);
    NumZeros = zeros(1,numel(BoundarySites));
    for b = 1:numel(BoundarySites)
        EdgeBonds = GraphArray{g}.Bonds(BoundarySites(b),:);
        NumZeros(b) = sum(EdgeBonds==0);
    end
    dZeros = NumZeros-EmptyBonds(g,:);
    Edge_cond = (sum(abs(dZeros)) == 0);
    if ~Edge_cond
        disp(['Number of empty bonds by corner: ' num2str(NumZeros)]);
        disp(['Expected zero bonds: ' num2str(EmptyBonds(g,:))]);
    end
    assert(Edge_cond,'Test failed: edge of lattice.');    
end

%% Test 2: Sublattice counting

for g = 1:numel(GraphArray)
    disp(['Testing Graph: ' GraphIDs{g}]);
    SLNum = max(GraphArray{g}.SLInds);
    Sublat_cond = (SLNum == SubLatNums(g));
    if ~Sublat_cond
        disp(['Maximum sublattice index: ' num2str(SLNum)]);
        disp(['Expected sublattice number: ' num2str(SubLatNums(g))]);
    end
    assert(Sublat_cond,'Test failed: number of sublattices.');    
end

%% Test 3: Number of invalid bonds

for g = 1:numel(GraphArray)
    disp(['Testing Graph: ' GraphIDs{g}]);
    ZeroBonds = sum(GraphArray{g}.Bonds(:)==0);
    Bonds_cond = (ZeroBonds == EmptyBondNums(g));
    if ~Bonds_cond
        disp(['Number of zeroed bonds: ' num2str(ZeroBonds)]);
        disp(['Expected number of empty bonds: ' num2str(EmptyBondNums(g))]);
    end
    assert(Bonds_cond,'Test failed: number of empty bonds')
end