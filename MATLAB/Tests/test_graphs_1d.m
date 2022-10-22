addpath(genpath("../Graphs"));

% Set some test cases with known outcomes.

% 1D: Nearest neighbour, next nearest neighbour, isolated, with both PBC
% and OBC.

G_1d_nn_pb = HypCub(10,1,1,1);
G_1d_nn_ob = HypCub(10,0,1,1);
G_1d_rn_pb = HypCub(10,1,-1,1);
G_1d_rn_ob = HypCub(10,0,-1,1);
G_1d_nnn_pb = HypCub(10,1,2,1);
G_1d_nnn_ob = HypCub(10,0,2,1);
G_1d_iso = HypCub(10,0,0,1);

GraphArray = {G_1d_nn_pb; G_1d_nn_ob; G_1d_rn_pb; G_1d_rn_ob; G_1d_nnn_pb; G_1d_nnn_ob; G_1d_iso};

GraphIDs = {'1D nearest neighbour, PBC'; '1D nearest neighbour, OBC';'1D reverse neighbour, PBC'; ...
    '1D reverse neighbour, OBC';'1D next nearest neighbour, PBC';'1D next nearest neighbour, OBC'; ...
    '1D isolated sites'};

BoundarySites = [1 10];
EmptyBonds = [0, 0;
              0, 1;
              0, 0;
              1, 0;
              0, 0;
              0, 1;
              0, 0];

SubLatNums = [1; 1; 1; 1; 2; 2; 10];

EmptyBondNums = [0; 1; 0; 1; 0; 2; 0];

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