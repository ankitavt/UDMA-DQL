clear all;
close all;
rng('shuffle');

playgroundfile = [];  % Run a new playground

if ~isempty(playgroundfile)
    load(playgroundfile);
    playgroundfile = 1;
    clear linkDQN;
    clear link;
end

% NMCruns = 1000; % Number of Mont Carlo runs
NMCruns = 10 % Number of Mont Carlo runs
Txpower = 40; % in dBm
% Normawgn = 4e-21; % AWGN per Hz
% No = -174; %5Noise PSD in dBm/Hz
% NF = 2.5; % Noise figure
% Bw = 180e3; % bandwidh for a resource block
Pn = -130; % Noise power in dBm
% PnArray = [-140:10:-50]; % Will analize performance as the thermal noise power changes
% PnArray = -50; % For debug purposes only
Pmax = 40;  % Maximum transmit power [dBm]
Pmin = -20; % Minimum transmit power [dBm]
PdBmChangeMin = 0.5; % If all allocated powers change by less than this, iteration ends

if isempty(playgroundfile) % Run next only if we haven't loaded a playground file (file name is empty)    
    ShadowStddB = 6; % Shadowing std in dB
    %ShadowStddB = 10; % Shadowing std in dB
    ShadowStd = db2pow(ShadowStddB); % Shadowing std in dB

    NeNB = 9; % Number of base stations, has to be a perfect square for a square playground
    Npue = 7; % Number of co-channel terminals in primary network (one in each eNB)
    %NCRTxRx = 40; % Number of cognitive radio transmitter-receiver pairs (can be placed anywhere in the playground)
    NCRTxRx = 2; % Number of cognitive radio transmitter-receiver pairs (can be placed anywhere in the playground)

    di = 200; % Distance between base stations
    CRMaxDist = 50; % Maximum distance between cognitive radio transmitter - receiver
end

BERtgt = 1e-3;
KforAdpMod = 1.5/-log(5*BERtgt);

PrimThrgRelChngLim = 0.05;
CompareRelChngLimAllCRs = PrimThrgRelChngLim * ones(NCRTxRx,1);

CRInitlSINR = 44*ones(1,NCRTxRx); % Initial target SINR for secondary network (in dB)
CRMinSINR = -16; % In [dB]
PmaxCR = 20;  % Maximum transmit power [dBm]
PminCR = -30; % Minimum transmit power [dBm]

% Available transmit powers, in effect the action setThis
ppossible = 10.^[-4:.25:-1];
ppossible = [0,ppossible]; % add no transmit action

policycnt=zeros(1,NMCruns);
pcounter = zeros (1,NMCruns);
% policycnt = 0;
% pcounter=0;
% % Table-based Q-learning parameters
NStates = 2;
dscnt = .2;
nExpPhases = 5;
%nExpPhases = 3;
ExpPhasesLen = 10;
ExpProb = 0.3;
%ExpProb = 0.6
ConvThr = 8;
%Tmprture = 5;
Tmprture = 2;
lrnrateUpdtExp=.51;
InertiaProb = 0.8;  % Probability of keeping policy (even when there may be other better policies)
TolrcLvl = 0.1;  % tolearance of Q value is relative to the maximum Q-value at the time

% Neural network Q-learning parameters
Nlayer = 4; % Number of layers
N_neuron_layer = [1,8,18,length(ppossible)]; % Neurons in each layer, one output neuron per action, input is state
LearnRatePol = 5e-5;
LearnRate = 0.03;
% LearnRate = 0.001;
mBatchSize = 5;
% EpLen = 150; % Learning epoch length
% nExpPhasesDQN = 650;
EpLen = 10; % Learning epoch length
nExpPhasesDQN = 5;
%ExpProbDQN = 0.1 % Probability of trying something different than the policy
ExpProbDQN = 0.3; % Probability of trying something different than the policy
TmprtureDQN = 3;
InertiaProbDQN = 0.25;  % Probability of keeping policy (even when there may be other better policies)
TolrcLvlDQN = 0;  % tolearance of Q value is relative to the maximum Q-value at the time
CParamDQN = 50; % How often to update target action Q values
%CParamDQN = 1; % How often to update target action Q values
useSoftMax = 0;
LearnRateDivBy = 1.1;
MaxLearnRateDivFac = 5000;
LenAllSimAnn = 5;
LenOneSimAnn = 1;
WrkSPcSvFile = ['network1.mat'];
diary sthtest;
diary on;
shutdownCR1Exp = 0;
shutdownCR2Exp = 0;
DiffLRate = 1;

if isempty(playgroundfile) % Run next only if we haven't loaded a playground file (file name is empty)    
    % Setup a grid of transmitters, each covering a square area (Fig. 1 in
    % paper). Code is a modified version of sinr_throughput_cellular_playground.m
    % All distances in meters
    % Distance resolution (main reason is to reduce error propagation ff_T_VehA_3retx_1ant
    % due to the presence of sqrt(3) and "x/3" in the distances
    dres = 0.1;

    % Measuring points resolution (or separation in meters)
    MPres = 5;

    % distance too close to eNodeB
    nearfld = 5;

    % Playground limits (in dots with resolution dres*dres)
    plygrnv = di * sqrt(NeNB)/dres;
    plygrnh = di * sqrt(NeNB)/dres;

    eNodeBlocH = ones(sqrt(NeNB),1) * [di/2/dres : di/dres : plygrnh];
    eNodeBlocV = eNodeBlocH';

    % I think that at this point dres will not be needed and we can do away with it
    eNodeBlocV = eNodeBlocV*dres;
    eNodeBlocH = eNodeBlocH*dres;

    % Define sample points on playground where to do measurements
    mpgridy = 0:MPres:plygrnv*dres;
    mpgridx = 0:MPres:plygrnh*dres;

    % Randomly choose eNBs that are active (those with a primary UR assigned)
    Act_eNB_indx = datasample(1:NeNB,Npue,'Replace',false)

    side_NeNB = sqrt(NeNB); % Number of eNBs in one side of the square playground
    Orig_eNB_1Dindx = reshape([1:NeNB],side_NeNB,side_NeNB);

    HalfSideLength = di * side_NeNB/2;
end

for ncn=1:NCRTxRx; PwrIndVec{ncn} = 1:length(ppossible)+1; end;
ppossiblePlusNoTx = [0 ppossible];
AllPwrComb = ppossiblePlusNoTx(combvec(PwrIndVec{:}) );

% Following performance results are from using lte_amc.m
T_AMC_PedB_3retx_1ant_BLER1em1 = [1e-10 0.0026    0.0088    0.0157    0.0228    0.0295    0.0367    0.0423    0.0458    0.0475    0.0869    0.1578    0.2582    ...
                                  0.3880    0.3960    0.6680    0.8602    0.8844    1.2429    1.2751    1.5460    1.8147    2.0947    2.1226    2.3688    ...
                                  2.4028    2.5582    2.6004    2.6268    2.6347    2.6374    2.6400];                              
ModType_AMC_PedB_3retx_1ant_BLER1em1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];

T_AMC_PedB_3retx_1ant_BLER1em2 = [1e-10 0.0026    0.0088    0.0157    0.0228    0.0295    0.0367    0.0423    0.0458    0.0475    0.0478    0.0480    0.0880    ...
                                  0.1600    0.2637    0.3992    0.4000    0.6859    0.6880    0.8960    1.2841    1.2880    1.5680    1.8480    2.1397    ...
                                  2.1440    2.1440    2.4296    2.4320    2.6347    2.6374    2.6400];
ModType_AMC_PedB_3retx_1ant_BLER1em2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];

SNR_PedB_3retx_1ant = [-40 -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45];

T_AMC_PedB_3retx_2x2_BLER1em1 = [1e-9 0.0024    0.0169    0.0325    0.0460    0.0619    0.0816    0.0934    0.0959    0.1598    0.2874    0.4938    0.7434    ...
                                 0.7520    1.0220    1.2837    1.6756    2.1054    2.1408    2.4309    2.9522    3.4760    4.0392    4.5570    4.8708    ...
                                 4.9989    5.0240    5.0240    5.0240    5.0240    5.0240    5.0240];
ModType_AMC_PedB_3retx_2x2_BLER1em1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
CQI_AMC_PedB_3retx_2x2_BLER1em1 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 9, 10, 11, 12, 13, 14, 15, 15, 15, 15, 15, 15, 15, 15];                             

T_AMC_PedB_3retx_2x2_BLER1em2 = [1e-9, 0.0024    0.0169    0.0325    0.0460    0.0619    0.0816    0.0934    0.0959    0.0960    0.1600    0.2880    0.4960    ...
                                 0.4960    0.7520    1.0240    1.2947    1.2960    1.6960    2.1440    2.4480    2.9760    3.5040    4.0759    4.0800    ...
                                 4.6240    5.0240    5.0240    5.0240    5.0240    5.0240    5.0240];
ModType_AMC_PedB_3retx_2x2_BLER1em2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];                             

%SNR_PedB_3retx_2x2 = [-16, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45];
SNR_PedB_3retx_2x2 = [-40, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 60];


%% Run the algorithm to calculate transmit power at the primary network
% Thrgpt = zeros(1,length(PnArray));
sinrUEFirspPass = [];
sinrPrimWithoutCR = [];
SINRPrimWithCR = [];
SINRPrimWithCRTXonModT0 = [];
ModTypeHist = [];
ModTypeHistNeighbr = [];
ThrgptRunNoCR = zeros(1,NMCruns);
RelThrChngRun=[];

for ncn=1:NStates; PwrIndCellRow{ncn} = 1:length(ppossible); end;
AllPolicies = ppossible(combvec(PwrIndCellRow{:}) );

% DonePhaseN = nExpPhases*ones(NCRTxRx,NMCruns); 

for u=1:NCRTxRx
    link(u).QmatTrackAppAvg=[];
    link(u).QmatTrackAppStd=[];
end

% SizeNextPolcySet = zeros(NMCruns * NCRTxRx,nExpPhasesDQN);

% ThrZero = 0;
mcrun = 1;
while mcrun <= NMCruns  % Main loop, a new network setup each time
    tic

%     for n = 1 : NCRTxRx
%         link(n).QPolOneExp = zeros(1,length(AllPolicies));
%     end
    
    if isempty(playgroundfile) % Run next only if we haven't loaded a playground file (file name is empty)    
        % ######## Start setup of networks (place nodes, calculate channel gains, etc.) #########
    
        % Randomly choose eNBs that are active (those with a primary UR assigned)
        Act_eNB_indx = datasample(1:NeNB,Npue,'Replace',false);

        % The whole next double for loop is basically to only get the area covered by active base stations
        % and get the channel gains for the covered points
        % Calculate SINR for each point of interest
        for cy = 1:length(mpgridy)
            for cx =  1:length(mpgridx)
                y = mpgridy(cy);
                x = mpgridy(cx);
                % Calculate distance to all eNodeBs
                DistX = abs(x-eNodeBlocH);
                DistY = abs(y-eNodeBlocV);
                % Process the distances to create wrap around effect and remove edges
                DistX(find(DistX > HalfSideLength)) = 2*HalfSideLength - DistX(find(DistX > HalfSideLength));
                DistY(find(DistY > HalfSideLength)) = 2*HalfSideLength - DistY(find(DistY > HalfSideLength));
                distances = sqrt(DistX.^2+DistY.^2);
                % Calculate received power from all eNodeBs
                % Txpower-30 converts dBm to dB, also convert distances to km
                % The paper H. W. Arnold, D. C. Cox and R. R. Murray, "Macroscopic diversity performance measured ...
                % in the 800-MHz portable radio communications environment," in IEEE Transactions on Antennas and Propagation, ... 
                % vol. 36, no. 2, pp. 277-281, Feb 1988, explains very well shadow loss
                GChannel = 10.^(.1*(-(128.1+37.6*log10(distances/1000))-10-ShadowStd*randn(size(distances))));
    %            GChannel = db2pow(-10*4.5*log10(distances/10) - ShadowStd*randn(size(distances)));
                RxPower = 10.^(.1*(Txpower-30)) .* GChannel;

                % Find index for eNodeB with max received power
                [clstsENBr,clstsENBc]=find(RxPower==max(max(RxPower)));

                % Option to max received power is eNodeB that is closest
                RxPowerVector=reshape(RxPower,1,numel(RxPower));
                [SortedRxPower,idxSortPrx]=sort(RxPowerVector,'descend');

                clstsENBr(2:end)=[]; % in case there is more than one answer
                clstsENBc(2:end)=[]; % in case there is more than one answer

                % Prepare channel gain info for saving
                GChannelVct = GChannel(clstsENBr,clstsENBc); % Place first the non-interfering channel
                GChannel(clstsENBr,clstsENBc) = -11; % Mark the location of the non-interfering channel
                GChannelVct = [GChannelVct GChannel(1:end)]; % Add the interfering channels as a vector

                ServENodeB(cy,cx) = idxSortPrx(1); % Record serving eNodeB

                if ~isempty(find(ServENodeB(cy,cx) == Act_eNB_indx)) % Check if this is a point covered by an active eNB (with no edges)
                    % Save channel gains (main and interfering ones)
                    GChannel_sv(cy,cx,1:length(GChannelVct)) = GChannelVct;
                end   
            end
        end

        % Get the actual location for each UE
        for n = 1 : Npue
            [CoveredAreaY,CoveredAreaX,values] = find(ServENodeB == Act_eNB_indx(n));
            rndI = randi([1 length(CoveredAreaY)],1,1);
            UElocationY(n) = CoveredAreaY(rndI); 
            UElocationX(n) = CoveredAreaX(rndI); 
            UElocationCoorX(n) = UElocationX(n) * MPres;  % AK Fixed bug 2/3/19 and on 5/30/19 added " * MPres"
            UElocationCoorY(n) = UElocationY(n) * MPres;  % AK Fixed bug 2/3/19 and on 5/30/19 added " * MPres"
       end

        % Deploy the coordinates for secondary, cognitive radio, network
        % Find first location for transmitters
        % AK Fixed bug 7/11/19: non-zero probability of colacating CR transmitter and primary receiver.
        % AK Fixed bug 7/11/19: non-zero probability of colacating CR transmitter and primary receiver.
        % AK Fixed bug 8/13/19: non-zero probability of colacating CR transmitter and primary transmitter.
        SameLoc = 1;
        while(SameLoc)
            SameLoc = 0;
            CRTxLocY = datasample(mpgridy,NCRTxRx,'Replace',false);
            CRTxLocX = datasample(mpgridx,NCRTxRx,'Replace',false);
            for CRTXcntr = 1:NCRTxRx
                % Check that a CR transmitter is not place at the same location of a primary receiver (it may happen with some probability)
                XCoorMatch = find(CRTxLocX(CRTXcntr)==UElocationCoorX);
                YCoorMatch = find(CRTxLocY(CRTXcntr)==UElocationCoorY);
                XCoorMatch(end+1:max([length(XCoorMatch),length(YCoorMatch)]))=NaN;
                YCoorMatch(end+1:max([length(XCoorMatch),length(YCoorMatch)]))=NaN;
                if(~isempty(XCoorMatch) && ~isempty(YCoorMatch) && ~isempty(intersect(XCoorMatch,YCoorMatch)))
                    SameLoc = 1;
                end
                % Check that a CR transmitter is not place at the same location of a primary transmitter either (it may happen with some probability)
                XCoorMatch = find(CRTxLocX(CRTXcntr)==eNodeBlocH(Act_eNB_indx));
                YCoorMatch = find(CRTxLocY(CRTXcntr)==eNodeBlocV(Act_eNB_indx));
                XCoorMatch(end+1:max([length(XCoorMatch),length(YCoorMatch)]))=NaN;
                YCoorMatch(end+1:max([length(XCoorMatch),length(YCoorMatch)]))=NaN;
                if(~isempty(XCoorMatch) && ~isempty(YCoorMatch) && ~isempty(intersect(XCoorMatch,YCoorMatch)))
                    SameLoc = 1;
                end
            end
        end

        % Find location for receivers within a range from each transmitter
        % and calculate some related channel gains
        for n = 1 : NCRTxRx
            % First see how close is the nearest other CR transmitter because
            % the receiver will need to be closer than this distance
            DistX = abs(CRTxLocX-CRTxLocX(n));
            DistY = abs(CRTxLocY-CRTxLocY(n));
            % Process the distances to create wrap around effect and remove edges
            DistX(find(DistX > HalfSideLength)) = 2*HalfSideLength - DistX(find(DistX > HalfSideLength));
            DistY(find(DistY > HalfSideLength)) = 2*HalfSideLength - DistY(find(DistY > HalfSideLength));
            distances = sqrt(DistX.^2+DistY.^2);
            distances(n) = CRMaxDist; % Since this distance has to be zero, I use it to also include the max. distance for receivers
            distCRfrTx = nearfld + (min(distances)-nearfld) * rand(1,NCRTxRx);
            AnglCRfrTx = 2 * pi * rand(1,length(NCRTxRx));    
            CRRxLocY = round(CRTxLocY + distCRfrTx * sin(AnglCRfrTx));
            CRRxLocX = round(CRTxLocX + distCRfrTx * cos(AnglCRfrTx));
           % Wrap coordinates around edges
            if CRRxLocY(n) > 2 * HalfSideLength;  CRRxLocY(n) = CRRxLocY(n) - 2 * HalfSideLength; end
            if CRRxLocY(n) < 0; CRRxLocY(n) = CRRxLocY(n) + 2 * HalfSideLength; end
            if CRRxLocX(n) > 2 * HalfSideLength;  CRRxLocX(n) = CRRxLocX(n) - 2 * HalfSideLength; end
            if CRRxLocX(n) < 0; CRRxLocX(n) = CRRxLocX(n) + 2 * HalfSideLength; end

            % Use loop to also calculate distance and channel gains from primary transmitters to CR
            % receiver and within CR network
            DistX = abs(eNodeBlocH(Act_eNB_indx)-CRRxLocX(n));
            DistY = abs(eNodeBlocV(Act_eNB_indx)-CRRxLocY(n));
            % Process the distances to create wrap around effect and remove edges
            DistX(find(DistX > HalfSideLength)) = 2*HalfSideLength - DistX(find(DistX > HalfSideLength));
            DistY(find(DistY > HalfSideLength)) = 2*HalfSideLength - DistY(find(DistY > HalfSideLength));
            distances = sqrt(DistX.^2+DistY.^2);
            % Calculate received power from all eNodeBs
            % Txpower-30 converts dBm to dB, also convert distances to km
            % The paper H. W. Arnold, D. C. Cox and R. R. Murray, "Macroscopic diversity performance measured ...
            % in the 800-MHz portable radio communications environment," in IEEE Transactions on Antennas and Propagation, ... 
            % vol. 36, no. 2, pp. 277-281, Feb 1988, explains very well shadow loss
            PrimtoCRChannel(n,1:Npue) = 10.^(.1*(-(128.1+37.6*log10(distances/1000))-10-ShadowStd*randn(size(distances))));
    %         PrimtoCRChannel(n,1:Npue) = db2pow(-10*4.5*log10(distances/10) - ShadowStd*randn(size(distances)));

            % Within CR next
            DistX = abs(CRTxLocX-CRRxLocX(n));
            DistY = abs(CRTxLocY-CRRxLocY(n));
            % Process the distances to create wrap around effect and remove edges
            DistX(find(DistX > HalfSideLength)) = 2*HalfSideLength - DistX(find(DistX > HalfSideLength));
            DistY(find(DistY > HalfSideLength)) = 2*HalfSideLength - DistY(find(DistY > HalfSideLength));
            distances = sqrt(DistX.^2+DistY.^2);
            % Calculate received power from all eNodeBs
            % Txpower-30 converts dBm to dB, also convert distances to km
            % The paper H. W. Arnold, D. C. Cox and R. R. Murray, "Macroscopic diversity performance measured ...
            % in the 800-MHz portable radio communications environment," in IEEE Transactions on Antennas and Propagation, ... 
            % vol. 36, no. 2, pp. 277-281, Feb 1988, explains very well shadow loss
            CRChannel(n,1:NCRTxRx) = 10.^(.1*(-(128.1+37.6*log10(distances/1000))-10-ShadowStd*randn(size(distances))));
    %        CRChannel(n,1:NCRTxRx) = db2pow(-10*4.5*log10(distances/10) - ShadowStd*randn(size(distances)));

            % Within CR transmitters next
            DistX = abs(CRTxLocX-CRTxLocX(n));
            DistY = abs(CRTxLocY-CRTxLocY(n));
            % Process the distances to create wrap around effect and remove edges
            DistX(find(DistX > HalfSideLength)) = 2*HalfSideLength - DistX(find(DistX > HalfSideLength));
            DistY(find(DistY > HalfSideLength)) = 2*HalfSideLength - DistY(find(DistY > HalfSideLength));
            distances = sqrt(DistX.^2+DistY.^2);
            % Calculate received power from all eNodeBs
            % Txpower-30 converts dBm to dB, also convert distances to km
            % The paper H. W. Arnold, D. C. Cox and R. R. Murray, "Macroscopic diversity performance measured ...
            % in the 800-MHz portable radio communications environment," in IEEE Transactions on Antennas and Propagation, ... 
            % vol. 36, no. 2, pp. 277-281, Feb 1988, explains very well shadow loss
            CRTXsChannel(n,1:NCRTxRx) = 10.^(.1*(-(128.1+37.6*log10(distances/1000))-10-ShadowStd*randn(size(distances))));
    %        CRChannel(n,1:NCRTxRx) = db2pow(-10*4.5*log10(distances/10) - ShadowStd*randn(size(distances)));

        end

        % Calculate distance and channel gains from CRs to primary receivers
        % and (to get primary transmitter receiver at CR transmitter with more power) to primary transmitters
        for n = 1:Npue
            DistX = abs(UElocationCoorX(n)-CRTxLocX);  % AK Fixed bug 2/3/19
            DistY = abs(UElocationCoorY(n)-CRTxLocY);  % AK Fixed bug 2/3/19
            % Process the distances to create wrap around effect and remove edges
            DistX(find(DistX > HalfSideLength)) = 2*HalfSideLength - DistX(find(DistX > HalfSideLength));
            DistY(find(DistY > HalfSideLength)) = 2*HalfSideLength - DistY(find(DistY > HalfSideLength));
            distances = sqrt(DistX.^2+DistY.^2);
            % Calculate received power from all eNodeBs
            % Txpower-30 converts dBm to dB, also convert distances to km
            % The paper H. W. Arnold, D. C. Cox and R. R. Murray, "Macroscopic diversity performance measured ...
            % in the 800-MHz portable radio communications environment," in IEEE Transactions on Antennas and Propagation, ... 
            % vol. 36, no. 2, pp. 277-281, Feb 1988, explains very well shadow loss
            CRtoPrimChannel(n,1:NCRTxRx) = 10.^(.1*(-(128.1+37.6*log10(distances/1000))-10-ShadowStd*randn(size(distances))));
    %         CRtoPrimChannel(n,1:NCRTxRx) = db2pow(-10*4.5*log10(distances/10) - ShadowStd*randn(size(distances)));
            DistX = abs(eNodeBlocH(Act_eNB_indx(n))-CRTxLocX);
            DistY = abs(eNodeBlocV(Act_eNB_indx(n))-CRTxLocY);
            % Process the distances to create wrap around effect and remove edges
            DistX(find(DistX > HalfSideLength)) = 2*HalfSideLength - DistX(find(DistX > HalfSideLength));
            DistY(find(DistY > HalfSideLength)) = 2*HalfSideLength - DistY(find(DistY > HalfSideLength));
            distances = sqrt(DistX.^2+DistY.^2);
            % Calculate received power from all eNodeBs
            % Txpower-30 converts dBm to dB, also convert distances to km
            % The paper H. W. Arnold, D. C. Cox and R. R. Murray, "Macroscopic diversity performance measured ...
            % in the 800-MHz portable radio communications environment," in IEEE Transactions on Antennas and Propagation, ... 
            % vol. 36, no. 2, pp. 277-281, Feb 1988, explains very well shadow loss
            PrimTxtoCRTxCh(n,1:NCRTxRx) = 10.^(.1*(-(128.1+37.6*log10(distances/1000))-10-ShadowStd*randn(size(distances))));
        end
        
        % ####### Finished setup of networks (place nodes, calculate channel gains, etc.) #######          

        PTx = Pmin + (Pmax-Pmin).*rand(1,length(Act_eNB_indx));  % Initial transmit power setting primary network
        PTxPrev = zeros(size(PTx));
        % PTx = 40*ones(size(PTx)); % For debug purposes only

        PTxCR = -Inf*ones(1,NCRTxRx);  % Initial transmit power [dBm] setting secondary network (equal to zero Watts because they are just listening) 

        % ####### Calculate transmit power at primary network #######
        % Update transmit power
        % Note that power update equation uses the channel gain and interference from same network. 
        AbsPChange = abs(PTx - PTxPrev); % To see power change from iteration to iteration (and decide when to stop)
        IterCnt = 1;
        while ( length(find(AbsPChange > PdBmChangeMin)) > 0  ) 
            [TotUEIntf,TotUERxTxP,sinrPrimWithoutCR,ThrgptPrim,IndxModTypeNoCR,CQI_No_CR] = Calc_PN_Metrics(side_NeNB,PTx,PTxCR,UElocationY,UElocationX,Pn,...
                                                GChannel,GChannel_sv,CRtoPrimChannel,Act_eNB_indx,SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,...
                                                ModType_AMC_PedB_3retx_2x2_BLER1em1,sinrPrimWithoutCR,CQI_AMC_PedB_3retx_2x2_BLER1em1);

            PTxPrev = PTx;
            % Calculate the updated primary power at one iteration
            PTx = UpdPrimPwrOneIter(Npue, Act_eNB_indx, UElocationX, UElocationY, GChannel_sv, TotUEIntf, Pn, Pmax, Pmin); 
            AbsPChange = abs(PTx - PTxPrev);
            IterCnt = IterCnt + 1;
        end
        % Calculate bunch of performance variables in primary network (without secondary network)
        [TotUEIntf,TotUERxTxP,sinrPrimWithoutCR,ThrgptPrimNoCR,IndxModTypeNoCR,CQI_No_CR] = Calc_PN_Metrics(side_NeNB,PTx,PTxCR,UElocationY,UElocationX,Pn,...
                                            GChannel,GChannel_sv,CRtoPrimChannel,Act_eNB_indx,SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,...
                                            ModType_AMC_PedB_3retx_2x2_BLER1em1,sinrPrimWithoutCR,CQI_AMC_PedB_3retx_2x2_BLER1em1);

        ThrgptRunNoCR(mcrun) = mean(ThrgptPrimNoCR); % Save result at this run

        % ----- Start for debugging purposes -----%
        NoiseAtPrmy =  db2pow(Pn-30); % Calculate noise at primary UEs
        sinrUE_No_CR = pow2db(  TotUERxTxP ./(TotUEIntf + NoiseAtPrmy) );
        PTx_No_CR = PTx;
        [ThrgptPrim_No_CR,IndxModType, CQI_No_CR] = GetThrghptandModType(sinrUE_No_CR,SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,ModType_AMC_PedB_3retx_2x2_BLER1em1,CQI_AMC_PedB_3retx_2x2_BLER1em1);
        % ----- End for debugging purposes -----%


        % New way of getting "closest" primary transmitter: Pick the one associated with largest channel gain
        for n = 1:NCRTxRx
            RxPoweratCRTx = PrimTxtoCRTxCh(:,n)'.*db2pow(PTx-30);
            ClosestCRTx(n) = min( find(RxPoweratCRTx == max(RxPoweratCRTx)) ); % Not actual eNB number but index into vector Act_eNB_indx of active eNB list
    %        ClosestCRTx(n) = min( find(PrimTxtoCRTxCh(:,n) == max(PrimTxtoCRTxCh(:,n))) ); % Not actual eNB number but index into vector Act_eNB_indx of active eNB list
        end

        ModTypeHist = [ModTypeHist IndxModTypeNoCR];
        ModTypeHistNeighbr = [ModTypeHistNeighbr IndxModTypeNoCR(unique((ClosestCRTx)))];

    end % end of skipping when loading playground file

    % ############# For debugging purposes ###############
    % ########## Plot Play Field ######
    % ################################################################    
%     PlotPlayField(CRTxLocX,CRTxLocY,CRRxLocX,CRRxLocY,eNodeBlocH,eNodeBlocV,UElocationCoorX,UElocationCoorY,Act_eNB_indx,NCRTxRx,Npue,mpgridx,mpgridy)
    % ######## End of for debugging purposes #############
    
    

    
    % ###### Add secondary network
         if (ThrgptRunNoCR(mcrun) ~= Inf) && ~isnan(ThrgptRunNoCR(mcrun))
            % Save workspace
%         clear all
%         load('network1.mat')
%         ExpProbDQN=0.009;
%         save(WrkSPcSvFile);
        clear linkDQN;
        clear link;        
        CRSINRTgt = CRInitlSINR;
        IndxModType = 4*ones(size(IndxModTypeNoCR));
        ModTypeDif = IndxModType - IndxModTypeNoCR;
        FirstPass = 1;

   
        
        % ###############################################
        % ############# Q-learning ######################
        % ###############################################

        % Q-learning initialization
        
        % Initialize table-based Q-learning
        SamePolcyCount = zeros(1,NCRTxRx);
        LinkDone = zeros(1,NCRTxRx);     
        for u=1:NCRTxRx
%            link(u).Qmat = ones(NStates,length(ppossible)); % Initialize Q-tables
            link(u).Qmat = rand(NStates,length(ppossible)); % Initialize Q-tables
            link(u).QmatEntry = zeros(NStates,length(ppossible)); % Initialize visit counter in Q-tables
            link(u).policy = ppossible(randi(length(ppossible),1,NStates)); 
            link(u).ActAvailbl = [ones(NStates,1) * ppossible , link(u).policy']; % Each row a different state
            link(u).PchooseActPMF = [ ExpProb/length(ppossible) : ExpProb/length(ppossible) : ExpProb , 1];
            link(u).RecentPast = zeros(NStates,ConvThr);
            link(u).QmatTrackAppNS = zeros(nExpPhases,ExpPhasesLen);
            link(u).QmatTrackApp = zeros(nExpPhases,ExpPhasesLen);
            link(u).Actions = ppossible;
            link(u).CountMinInstnc = zeros(size(link(u).Qmat));
            link(u).oldpolicy = ppossible(randi(length(ppossible),1,NStates));
            link(u).PolicySv = [];
        end
        xt=ones(1,NCRTxRx); % Initial state, with no CR transmission yet, the primary throughput can't have changed

        RelThrChngPhase = zeros(1,nExpPhases);
        p = zeros(nExpPhases,NCRTxRx);
        
        % ---------------------------------------------
        % ---------------------------------------------
        % +++++++++  Table-based Q-learning  ++++++++++
        % ---------------------------------------------
        % ---------------------------------------------
        c=1;
       % choose action
            % choose action (and execute action)
            for t = 1:ExpPhasesLen
                rndraw=rand(1,NCRTxRx);
                for u = 1:NCRTxRx
                    if ~LinkDone(u)


                        % Choose actions
                        pisuser = link(u).PchooseActPMF;  % pisuser is a vector with a "CDF" for actions
                        actindx = max(find(rndraw(u)>pisuser))+1; % pick action index from the random draw in rndraw
                        if(isempty(actindx))
                            actindx = 1;
                        end
                        actindxsv(u) = find(link(u).Actions == link(u).ActAvailbl(xt(u),actindx)); 
                        p(c,u) = link(u).ActAvailbl(xt(u),actindxsv(u)); % Action for link u at interation c
                    else
                        p(c,u) = link(u).policy(xt(u));
                    end  

                end    
                PTxRes = p(c,:); % Power at CRs
                PTxCR = pow2db(PTxRes) + 30; % Power at CRs converted to dBm 

%                CRtoPrimChannel = ones(size(CRtoPrimChannel)); % For debugging: all CR should NOT transmit with this channel

%                 % Get new state (indicator function of exceeding relative throughput change in PN
                NoiseAtPrmy = db2pow(PTxCR-30) * CRtoPrimChannel' + db2pow(Pn-30); % Calculate noise at primary UEs, secondary networks appears as noise
                sinrUE = pow2db(  TotUERxTxP ./(TotUEIntf + NoiseAtPrmy) );
                % Get throughput and modulation index in primary network
                [ThrgptPrim,IndxModType] = GetThrghptandModType(sinrUE,SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,ModType_AMC_PedB_3retx_2x2_BLER1em1,CQI_AMC_PedB_3retx_2x2_BLER1em1);
%                 find( ModTypeDif(ClosestCRTx) < 0 )) = CRSINRTgt(find( ModTypeDif(ClosestCRTx) < 0 )
                RelThrChng = abs(ThrgptPrimNoCR - ThrgptPrim)./ ThrgptPrimNoCR;
                % Receive new state
                xtp1 = ones(1,NCRTxRx);
                xtp1(RelThrChng(ClosestCRTx) > PrimThrgRelChngLim) = 2;
%                 if(max(RelThrChng) > PrimThrgRelChngLim)
%                     xtp1 = xtp1 * 2;
%                 end    
% 
                % compute reward (function of CR throughput)
                for n = 1:NCRTxRx
                    IntfCh = CRChannel(n,:);
                    IntfCh(n) = [];
                    IntfPw = db2pow(PTxCR-30);
                    IntfPw(n) = [];
                    SINRDen(n) = IntfCh * IntfPw' + PrimtoCRChannel(n,:) * db2pow(PTx-30)' + db2pow(Pn-30); % Primary network appears as noise
                end
                SINRCRRx = diag(CRChannel)' .* db2pow(PTxCR-30) ./ SINRDen;
                % Calculate throughput at CRs
                ThrgptCR = GetThrghptandModType(pow2db(SINRCRRx),SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,ModType_AMC_PedB_3retx_2x2_BLER1em1,CQI_AMC_PedB_3retx_2x2_BLER1em1);
                ThrgptCR(PTxRes == 0) = 0;
%                reward = ThrgptCR;
%                reward(xtp1 == 2) = -1;
                reward = 2.^ThrgptCR;
                reward(xtp1 == 2) = 0;
                rewardsv(:,c*ExpPhasesLen+t) = reward;
                % Update Q-tables
                for u = 1:NCRTxRx
                    if ~LinkDone(u)
                        % Count visits to Q-table entry
                        link(u).QmatEntry(xt(u),actindxsv(u)) = link(u).QmatEntry(xt(u),actindxsv(u)) + 1;
                        lrnrate = 1/link(u).QmatEntry(xt(u),actindxsv(u))^lrnrateUpdtExp;

                        % Update Q-table
                        DeltaQ = lrnrate * ( reward(u) + dscnt*max(link(u).Qmat(xtp1(u),:)) - link(u).Qmat(xt(u),actindxsv(u)) );
                        link(u).Qmat(xt(u),actindxsv(u)) = link(u).Qmat(xt(u),actindxsv(u)) + DeltaQ;
                        % Track update to Q-value as approximation progresses (measuring convergence of approximation)
%                         PrevNSmpl = link(u).QmatTrackAppNS(c,link(u).QmatEntry(xt(u),actindxsv(u)));
%                         PrevCnt = link(u).QmatTrackApp(c,link(u).QmatEntry(xt(u),actindxsv(u)));
%                         link(u).QmatTrackAppNS(c,link(u).QmatEntry(xt(u),actindxsv(u))) = PrevNSmpl + 1;
%                         link(u).QmatTrackApp(c,link(u).QmatEntry(xt(u),actindxsv(u))) = (PrevNSmpl*PrevCnt + abs(DeltaQ))/(PrevNSmpl + 1) ;
                    end
                end    
                xt=xtp1;
            end % End of one exploration phase, loop runs ExpPhasesLen times
% 
% %             figure();plotyy(1:length(ppossible),link(u).Qmat',1:length(ppossible),UndLayCond);grid  % For debug only
% 
            for u=1:NCRTxRx
                % Find set of valid next policies
                for xcnt = 1:NStates
                    % AK: This part next to check
                    maxQEntrySt = max(link(u).Qmat(xcnt,link(u).QmatEntry(xcnt,:) > 2)); % Don't consider Q-values approximated from few samples
                    if isempty(maxQEntrySt)
                        ValActInd{xcnt} = 1;
                        maxQEntrySt = link(u).Qmat(xcnt,1);
                    else    
                        ValActInd{xcnt} = find( link(u).Qmat(xcnt,:) == maxQEntrySt);
                    end
%                    ValActInd{xcnt} = find( link(u).Qmat(xcnt,:) >= max(link(u).Qmat(xcnt,:)) - TolrcLvl);
                    ValActInd{xcnt} = find( link(u).Qmat(xcnt,:) >= max(link(u).Qmat(xcnt,:)) * (1 - TolrcLvl) );
                    % ====== For debug =====
%                     link(u).indxPpossiblePol(xcnt) = find(link(u).Actions == link(u).policy(xcnt));
%                     link(u).QmatPol(xcnt,c) = link(u).Qmat(xcnt,link(u).indxPpossiblePol(xcnt));
                    link(u).MaxQmat(xcnt,c) = maxQEntrySt;  % For Debug
                    link(u).MaxQmatEnt(xcnt,c) = mean( link(u).QmatEntry(xcnt,find( link(u).Qmat(xcnt,:) == maxQEntrySt)) );
%                     link(u).PolQmatEnt(xcnt,c) = link(u).QmatEntry(xcnt,link(u).indxPpossiblePol(xcnt));
                    link(u).Qmatsv(xcnt+(c-1)*NStates,:) = link(u).Qmat(xcnt,:); 
                    link(u).QmatEntrysv(xcnt+(c-1)*NStates,:) = link(u).QmatEntry(xcnt,:); 
                    % ====== End Debug ======
                end
% 
                NextPolcySet = link(u).Actions(combvec(ValActInd{:}) ); % Each column is a policy (as transmit power level), rows are states 
                if size(NextPolcySet,1) == 1; NextPolcySet = NextPolcySet'; end; 

                % See if current policy belongs to valid next policies set
                IndxCurPolInPlcySet = find(sum(bsxfun(@eq,NextPolcySet,link(u).policy'))==length(link(u).policy')); % Current policy is NextPolcySet(:,IndxCurPolInPlcySet)
                if ~isempty(IndxCurPolInPlcySet) 
                    % Belongs
                    link(u).policy = link(u).policy; % Next policy is current policy (no changes)
                else
                    % Current policy is not as good as those in NextPolcySet
                    ChoosePolcyPMF = [ (1-InertiaProb)/size(NextPolcySet,2) : (1-InertiaProb)/size(NextPolcySet,2) : 1-InertiaProb , 1];
                    nextpolindx = max(find(rand(1) > ChoosePolcyPMF))+1; % pick action index from the random draw in rndraw
                    if(isempty(nextpolindx))
                        nextpolindx = 1;
                    end
                    if nextpolindx == length(ChoosePolcyPMF)
                        link(u).policy = link(u).policy; % Next policy is current policy (no changes)
                    else
                        link(u).policy = NextPolcySet(:,nextpolindx)'; 
                    end    
                end  % end of checking if current policy belongs to valid next policies set

                link(u).ActAvailbl = [ones(NStates,1) * ppossible , link(u).policy']; % Each row a different state, AK bug fix 10/18/2019
            end
% 
            for n = 1:NCRTxRx
                PTxRes(n) = link(n).policy(1); % Power at CRs
%                 link(n).PolicyIndx = find(ismember(AllPolicies',link(n).policy,'rows')==1);
                link(n).PolicySv(c,:) = link(n).policy;
            end    
            PTxCR = pow2db(PTxRes) + 30; % Power at CRs converted to dBm 
            NoiseAtPrmy = PTxRes * CRtoPrimChannel' + db2pow(Pn-30); % Calculate noise at primary UEs, secondary networks appears as noise
            sinrUE = pow2db(  TotUERxTxP ./(TotUEIntf + NoiseAtPrmy) );
            % Get throughput and modulation index in primary network
            [ThrgptPrim,IndxModType] = GetThrghptandModType(sinrUE,SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,ModType_AMC_PedB_3retx_2x2_BLER1em1,CQI_AMC_PedB_3retx_2x2_BLER1em1);
            for n = 1:NCRTxRx
                IntfCh = CRChannel(n,:);
                IntfCh(n) = [];
                IntfPw = db2pow(PTxCR-30);
                IntfPw(n) = [];
                SINRDen(n) = IntfCh * IntfPw' + PrimtoCRChannel(n,:) * db2pow(PTx-30)' + db2pow(Pn-30); % Primary network appears as noise
            end
            SINRCRRx = diag(CRChannel)' .* db2pow(PTxCR-30) ./ SINRDen;
            % Calculate throughput at CRs
            ThrgptCRStd = GetThrghptandModType(pow2db(SINRCRRx),SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,ModType_AMC_PedB_3retx_2x2_BLER1em1,CQI_AMC_PedB_3retx_2x2_BLER1em1);

            RelThrChng = abs(ThrgptPrimNoCR - ThrgptPrim)./ ThrgptPrimNoCR;
            % record worst case of relative throughput case
            RelThrChngPhase(c) = max(RelThrChng);   
            CRThrghFromPolcy(c,:) = ThrgptCRStd;
% 
            c = c + 1
        end  % End of all exploration phases

            
        % -------------------------------------------------------
        % -------------------------------------------------------
        % ++++++++++  End of Table-based Q-learning  ++++++++++++
        % -------------------------------------------------------
        % -------------------------------------------------------
             
        % ---------------------------------------------
        % ---------------------------------------------
        % +++++++++  Deep Q-learning Network  +++++++++
        % ---------------------------------------------
        % ---------------------------------------------  

        % Initialize neural network for Deep Q-Learning
        % Available activation functions (in order as entered):
        % 1)ReLU, 2)linear, 3)saturated linear, 4)symmetric saturating linear,
        % 5)Log sigmoid 6)Hyperbolic tangent sigmoid
        Act_Fun = {@poslin,@purelin,@satlin,@satlins,@logsig,@tansig};
        Act_Fun_Drvtv = {@(x)x>0 ; @(x)ones(size(x)) ; @(x)ones(size(x))-((x<0)|(x>=1)).*ones(size(x)) ; ...
                         @(x)ones(size(x))-((x<-1)|(x>=1)).*ones(size(x)) ; ... 
                         @(x)logsig(x).*(1-logsig(x)) ; @(x)tansig(x).^2};
        for u = 1 : NCRTxRx
            linkDQN(u).my_nnEv.Nlayer = Nlayer;
            linkDQN(u).my_nnEv.Act_Fun = Act_Fun;
            linkDQN(u).my_nnEv.Act_Fun_Drvtv = Act_Fun_Drvtv;
            % Define activation functions
            linkDQN(u).my_nnEv.actFunIndx{1} = [];
            linkDQN(u).my_nnEv.actFunIndx{2} = 3;
            linkDQN(u).my_nnEv.actFunIndx{3} = 3;
            linkDQN(u).my_nnEv.actFunIndx{4} = 1;


            % Initialize weights matrices,  activation and bias vectors for DQN
            linkDQN(u).my_nnEv.a{1} = zeros( N_neuron_layer(1) , 1 );  
            for lay = 2:Nlayer
                % One matrix for each layer after first one
                % 'lay' index in linkDQN(u).my_nnEv.w{lay} is for layer "to the right". Same for biases.
                % Row is "TO" and column is "FROM" connection
                linkDQN(u).my_nnEv.w{lay} = rand( N_neuron_layer(lay) , N_neuron_layer(lay-1) );  % Weights
                linkDQN(u).my_nnEv.b{lay} = rand( N_neuron_layer(lay) , 1);  % Biases
                linkDQN(u).my_nnEv.a{lay} = zeros( N_neuron_layer(lay) , 1 );  % Outputs (after activation function)
            end        
            linkDQN(u).Qmat = zeros(NStates,length(ppossible)); % Initialize Q-tables
            linkDQN(u).QmatEntry = zeros(NStates,length(ppossible)); % Initialize visit counter in Q-tables
%           link(u).ActAvailbl = [ones(NStates,1) * ppossible , link(u).policy']; % Each row a different state
            linkDQN(u).PchooseActPMF = zeros(size(ppossible));

            % Initialization of some research variables

            linkDQN(u).SaveRwrdIndx(1:length(ppossible)) = 1; % Index (pointer) where to save reward
            linkDQN(u).SaveRwrd(1:length(ppossible),1) = 0; % Initialize storage of reward
            linkDQN(u).QmatSv = zeros(nExpPhasesDQN*NStates*EpLen,length(ppossible));  % Save the Q-values
            linkDQN(u).QmatSvPol = zeros(nExpPhasesDQN*EpLen,NStates);  % Save the Q-values of the policy    

        end  
        
        xt=ones(1,NCRTxRx); % Initial state, with no CR transmission yet, the primary throughput can't have changed
        xtsv = xt(1);

        actionDQNindexSv=[];

        for u=1:NCRTxRx
            % Make forward pass to get vector of Q-values for each state. 
            for  xcnt = 1:NStates
                linkDQN(u).my_nnEv = make_nn_fwd_pass(xcnt,linkDQN(u).my_nnEv);
                linkDQN(u).Qmat(xcnt,:) = linkDQN(u).my_nnEv.a{Nlayer};
            end
            linkDQN(u).Actions = ppossible;
            linkDQN(u).QmatTrgt = linkDQN(u).Qmat;  % Update target Q-values
            linkDQN(u).policy = ppossible(randi(length(ppossible),1,NStates)); 
%             linkDQN(u).PolicySv(1,1:NStates) = linkDQN(u).policy;
            [RowInd,PolIndxVal]=find(bsxfun(@eq,[ppossible;ppossible],linkDQN(u).policy'));  % AK 10/21/2019
            linkDQN(u).policyAsIndx(RowInd) = PolIndxVal; % AK 10/21/2019
            linkDQN(u).PolicySv(1,1:NStates) = linkDQN(u).policyAsIndx;  % AK 10/21/2019
            linkDQN(u).ActAvailbl = [ones(NStates,1) * ppossible , linkDQN(u).policy']; % Each row a different state
            linkDQN(u).LearnRateVec = LearnRate * ones(length(ppossible),1);
            linkDQN(u).PchooseActPMF = [ ExpProbDQN/length(ppossible) : ExpProbDQN/length(ppossible) : ExpProbDQN , 1];

            % Next for simulated annealing component
            for lay = 2:Nlayer
                % One matrix for each layer after first one
                % 'lay' index in linkDQN(u).my_nnEv.w{lay} is for layer "to the right". Same for biases.
                % Row is "TO" and column is "FROM" connection
                linkDQN(u).my_nnEv.wSv{lay} = linkDQN(u).my_nnEv.w{lay};  % Weights
                linkDQN(u).my_nnEv.bSv{lay} = linkDQN(u).my_nnEv.b{lay};  % Biases
            end    
            for  xcnt = 1:NStates            
%                 linkDQN(u).QSimAnnPolBef(xcnt) =  linkDQN(u).Qmat(xcnt,linkDQN(u).policyAsIndx(xcnt));
                linkDQN(u).QSimAnnPolBef(xcnt) =  -Inf;
            end    
            linkDQN(u).QSimAnnPolAsInd = linkDQN(u).policyAsIndx;
            linkDQN(u).QSimAnnPolicy = linkDQN(u).policy;
        end

        c=1;
        while c <= nExpPhasesDQN
            actionDQNindex = zeros(EpLen*mBatchSize,NCRTxRx);
            actionDQN = zeros(EpLen*mBatchSize,NCRTxRx);

            [mcrun  c]

            for ep = 1:EpLen  % A learning episode is also an exploration phase

                for u=1:NCRTxRx
                    % Initialize accumulators for minibatches weights/biases update
                    for lay = 2:Nlayer
                        linkDQN(u).AccW{lay} = zeros( N_neuron_layer(lay) , N_neuron_layer(lay-1) );  
                        linkDQN(u).AccB{lay} = zeros( N_neuron_layer(lay) , 1);  
                    end
                end
                for mm = 1:mBatchSize
                    rndraw=rand(1,NCRTxRx);
                    for u=1:NCRTxRx
                        if useSoftMax
                            linkDQN(u).PchooseActPMF(end) = [];
                            linkDQN(u).PchooseAct = exp(linkDQN(u).Qmat(xt(u),:)/TmprtureDQN); % Oldest way of doing this, and maybe incorrect
                            linkDQN(u).PchooseAct = linkDQN(u).PchooseAct/sum(linkDQN(u).PchooseAct);
                            linkDQN(u).PchooseActPMF(1) = linkDQN(u).PchooseAct(1); 
                            for cnt = 2:length(ppossible)
                                linkDQN(u).PchooseActPMF(cnt) = linkDQN(u).PchooseActPMF(cnt-1) + linkDQN(u).PchooseAct(cnt);
                            end       
                            linkDQN(u).PchooseActPMF = [linkDQN(u).PchooseActPMF * ExpProbDQN , 1];
                        end
                        % Choose actions
                        pisuser = linkDQN(u).PchooseActPMF;  % pisuser is a vector with a "CDF" for actions
                        actindx = max(find(rndraw(u)>pisuser))+1; % pick action index from the random draw in rndraw
                        if(isempty(actindx))
                            actindx = 1;
                        end

                        actionDQNindex((ep-1)*mBatchSize+mm,u) = min(find(ppossible == linkDQN(u).ActAvailbl(xt(u),actindx))); % Action index for link u at minibatch index mm

                        actionDQN((ep-1)*mBatchSize+mm,u) = linkDQN(u).ActAvailbl(xt(u),actionDQNindex((ep-1)*mBatchSize+mm,u)); % Action for link u at minibatch index mm
                    end

                    PTxRes = actionDQN((ep-1)*mBatchSize+mm,:); % Power at CRs
                    PTxCR = pow2db(PTxRes) + 30; % Power at CRs converted to dBm 

                    NoiseAtPrmy = db2pow(PTxCR-30) * CRtoPrimChannel' + db2pow(Pn-30); % Calculate noise at primary UEs, secondary networks appears as noise
                    sinrUE = pow2db(  TotUERxTxP ./(TotUEIntf + NoiseAtPrmy) );
                    % Get throughput and modulation index in primary network
                    [ThrgptPrim,IndxModType] = GetThrghptandModType(sinrUE,SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,ModType_AMC_PedB_3retx_2x2_BLER1em1,CQI_AMC_PedB_3retx_2x2_BLER1em1);
                    RelThrChng = abs(ThrgptPrimNoCR - ThrgptPrim)./ ThrgptPrimNoCR;
                    % Receive new state
                    xtp1 = ones(1,NCRTxRx);
%                     xtp1(RelThrChng(ClosestCRTx) > PrimThrgRelChngLim) = 2;
                    if( ~isempty(find(RelThrChng(ClosestCRTx) > PrimThrgRelChngLim)) )
                        xtp1 = xtp1 * 2;
                    end    
    %                 if(max(RelThrChng) > PrimThrgRelChngLim)
    %                     xtp1 = xtp1 * 2;
    %                 end    

                    % compute reward (function of CR throughput)
                    for n = 1:NCRTxRx
                        IntfCh = CRChannel(n,:);
                        IntfCh(n) = [];
                        IntfPw = db2pow(PTxCR-30);
                        IntfPw(n) = [];
                        SINRDen(n) = IntfCh * IntfPw' + PrimtoCRChannel(n,:) * db2pow(PTx-30)' + db2pow(Pn-30); % Primary network appears as noise
                    end
                    SINRCRRx = diag(CRChannel)' .* db2pow(PTxCR-30) ./ SINRDen;
                    % Calculate throughput at CRs
                    ThrgptCR = GetThrghptandModType(pow2db(SINRCRRx),SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,ModType_AMC_PedB_3retx_2x2_BLER1em1,CQI_AMC_PedB_3retx_2x2_BLER1em1);
                    ThrgptCR(PTxRes == 0) = 0;
                    reward = 1e6.^(sum(ThrgptCR)*ones(size(ThrgptCR))); % Most recent
                    reward(xtp1 == 2) = 0;

                    for u=1:NCRTxRx
                        % Calculate target
%                         linkDQN(u).trgt = linkDQN(u).QmatTrgt(xt(u),:)';
                        linkDQN(u).trgt = linkDQN(u).Qmat(xt(u),:)';
%                         linkDQN(u).trgt( actionDQNindex((ep-1)*mBatchSize+mm,u) ) = reward(u) + dscnt*max(linkDQN(u).Qmat(xtp1(u),:));
                        linkDQN(u).trgt( actionDQNindex((ep-1)*mBatchSize+mm,u) ) = reward(u) + dscnt*max(linkDQN(u).QmatTrgt(xtp1(u),:));
                        linkDQN(u).QmatEntry(xt(u),actionDQNindex((ep-1)*mBatchSize+mm,u)) = linkDQN(u).QmatEntry(xt(u),actionDQNindex((ep-1)*mBatchSize+mm,u)) + 1; % For debug

                        % Forward pass
                        linkDQN(u).my_nnEv = make_nn_fwd_pass(xt(u),linkDQN(u).my_nnEv);

                        % Backpopagation (assuming mean squared error)
                        [linkDQN(u).AccW,linkDQN(u).AccB] = bckprop_mse_updtw(xt(u),linkDQN(u).trgt,linkDQN(u).my_nnEv,linkDQN(u).AccW,linkDQN(u).AccB);

                        % Track rewards each time to research variability of environment
                        linkDQN(u).SaveRwrd(actionDQNindex((ep-1)*mBatchSize+mm,u),linkDQN(u).SaveRwrdIndx(actionDQNindex((ep-1)*mBatchSize+mm,u))) = reward(u);
                        linkDQN(u).SaveRwrdIndx(actionDQNindex((ep-1)*mBatchSize+mm,u)) = linkDQN(u).SaveRwrdIndx(actionDQNindex((ep-1)*mBatchSize+mm,u)) + 1;

%                       if ( linkDQN(u).ActAvailbl(xt(u),actindx)==linkDQN(u).policyAsIndx)
%                         paverage= zeros(1,NMCruns);
                        pestimate = zeros(1,NMCruns);
                        if nExpPhasesDQN >=15
                            if  linkDQN(u).ActAvailbl(xt(u),actindx) == linkDQN(u).ActAvailbl(xt(u),actionDQNindex((ep-1)*mBatchSize+mm,u))
                                 (mcrun) == policycnt (mcrun)+1;

                                if reward(u) == 0
                                pcounter (mcrun) = pcounter (mcrun)+1;              
                                end
                                pestimate = pcounter./policycnt;
                                paverage = pestimate(pestimate~=0);
                            end
                        end



                    end

                    xt=xtp1;
                    xtsv(end+1) = xt(1);


                end  % End mini-batch size

                for u=1:NCRTxRx

                    % Update weights and biases, gradient descent
                    for lay = 1:Nlayer
                        LearnRateArrayW = LearnRate * ones(size( linkDQN(u).AccW{lay}));
                        LearnRateArrayB = LearnRate * ones(size( linkDQN(u).AccB{lay}));
                        if lay == Nlayer
                            LearnRateArrayW = repmat(linkDQN(u).LearnRateVec, 1, size(linkDQN(u).AccW{lay},2) );
                            LearnRateArrayB = linkDQN(u).LearnRateVec;
                        end
                        linkDQN(u).my_nnEv.w{lay} = linkDQN(u).my_nnEv.w{lay} - LearnRateArrayW / mBatchSize .* linkDQN(u).AccW{lay};
                        linkDQN(u).my_nnEv.b{lay} = linkDQN(u).my_nnEv.b{lay} - LearnRateArrayB / mBatchSize .* linkDQN(u).AccB{lay};
                    end

                    % Make forward pass to get vector of Q-values for each state.policycnt
                    for  xcnt = 1:NStates
                        linkDQN(u).my_nnEv = make_nn_fwd_pass(xcnt,linkDQN(u).my_nnEv);
                        linkDQN(u).Qmat(xcnt,:) = linkDQN(u).my_nnEv.a{Nlayer};                    
                        linkDQN(u).QmatSv(xcnt+(ep-1)*NStates+(c-1)*NStates*EpLen,:) = linkDQN(u).Qmat(xcnt,:);  % Save the Q-values
                        linkDQN(u).QmatSvPol(ep+(c-1)*EpLen,xcnt) = linkDQN(u).Qmat(xcnt,linkDQN(u).PolicySv(c,xcnt));  % Save the Q-values of the policy
                    end
                    if mod(ep+(c-1)*EpLen,CParamDQN) == 0
                        linkDQN(u).QmatTrgt = linkDQN(u).Qmat;  % Update target Q-values
                    end

                    % Track error
                    linkDQN(u).Relerror_epoch(mcrun,ep+(c-1)*EpLen) = norm(linkDQN(u).trgt-linkDQN(u).my_nnEv.a{Nlayer})^2/norm(linkDQN(1).trgt)^2;
                end    
            end % End learning episode (also an exploration phase), end of one exploration phase, loop runs ExpPhasesLen times


            for u=1:NCRTxRx
                % Find set of valid next policies
%                linkDQN(u).QmatSTD = movstd(linkDQN(u).QmatSv(1:2:end,:),EpLen); % Calculate standard deviation of Q-values (think of having AWGN error from estimate
                linkDQN(u).QmatSTD = movstd(linkDQN(u).QmatSv(1+(c-1)*NStates*EpLen:2:2+(EpLen-1)*NStates+(c-1)*NStates*EpLen,:),EpLen); % Calculate standard deviation of Q-values (think of having AWGN error from estimate
                maxSTDQ = max(linkDQN(u).QmatSTD(end,:)); % Get largest standard deviation
%                TolrcLvlDQN = 3*maxSTDQ; % Tolerance is 3 times sigma (standard deviation) for maximum sigma
                for xcnt = 1:NStates
                    TolrcLvlDQN = 3*linkDQN(u).QmatSTD(end,find( linkDQN(u).Qmat(xcnt,:) == max(linkDQN(u).Qmat(xcnt,:)))); % Tolerance is 3 times sigma (standard deviation) for maximum sigma
                    ValActInd{xcnt} = find( linkDQN(u).Qmat(xcnt,:) >= max(linkDQN(u).Qmat(xcnt,:)) - TolrcLvlDQN);
                    linkDQN(u).TolrcBorder(xcnt,c) = max(linkDQN(u).Qmat(xcnt,:)) - TolrcLvlDQN;
                    if linkDQN(u).TolrcBorder(xcnt,c) < 0
                        linkDQN(u).TolrcBorder(xcnt,c)
                    end    
%                    ValActInd{xcnt} = find( linkDQN(u).Qmat(xcnt,:) >= max(linkDQN(u).Qmat(xcnt,:)) * (1 - TolrcLvlDQN) );
                end

                NextPolcySet = ppossible(combvec(ValActInd{:}) ); % Each column is a policy (as transmit power level), rows are states 
                if size(NextPolcySet,1) == 1; NextPolcySet = NextPolcySet'; end; 
                linkDQN(u).SizeNextPolcySetSv(c) = size(NextPolcySet,1);

%                     % Current policy is not as good as those in NextPolcySet
                ChoosePolcyPMF = [ (1-InertiaProbDQN)/size(NextPolcySet,2) : (1-InertiaProbDQN)/size(NextPolcySet,2) : 1-InertiaProbDQN , 1];
                nextpolindx = max(find(rand(1) > ChoosePolcyPMF))+1; % pick action index from the random draw in rndraw
                if(isempty(nextpolindx))
                    nextpolindx = 1;
                end
                if nextpolindx == length(ChoosePolcyPMF)
                    linkDQN(u).policy = linkDQN(u).policy; % Next policy is current policy (no changes)
                else
                    linkDQN(u).policy = NextPolcySet(:,nextpolindx)'; 
                end    
%                 end  % end of checking if current policy belongs to valid next policies set
                linkDQN(u).ActAvailbl = [ones(NStates,1) * ppossible , linkDQN(u).policy'];
                [RowInd,PolIndxVal]=find(bsxfun(@eq,[ppossible;ppossible],linkDQN(u).policy'));  % AK 10/21/2019
                linkDQN(u).policyAsIndx(RowInd) = PolIndxVal; % AK 10/21/2019
                linkDQN(u).PolicySv(c+1,1:NStates) = linkDQN(u).policyAsIndx;  % AK 10/21/2019
%                 linkDQN(u).PsdoSNRQPol = mean((10*log10(linkDQN(u).QmatSv(EpLen:2:end,:).^2./movstd(linkDQN(u).QmatSv(EpLen:2:end,:),EpLen).^2)));
%                 linkDQN(u).PolicySv(c+1,1:NStates) = linkDQN(u).policy;
                linkDQN(u).LearnRateVec(linkDQN(u).policyAsIndx(1),:) = linkDQN(u).LearnRateVec(linkDQN(u).policyAsIndx(1),:) / LearnRateDivBy;
                linkDQN(u).LearnRateVec(linkDQN(u).LearnRateVec < LearnRate/MaxLearnRateDivFac) = LearnRate/MaxLearnRateDivFac;
            end % end loop going through all CRs


            [linkDQN(1).policyAsIndx(1),linkDQN(2).policyAsIndx(1),BestCRPwrIndxSv(mcrun,:),[linkDQN(1).policyAsIndx(1),linkDQN(2).policyAsIndx(1)]-BestCRPwrIndxSv(mcrun,:),BestThrgptEachCR(mcrun,:)]

            % Do "simulated annealing" work
            if (c <= LenAllSimAnn) && (mod(c,LenOneSimAnn) == 0)
                for u=1:NCRTxRx
                    if linkDQN(u).Qmat(1, linkDQN(u).policyAsIndx(1)) > linkDQN(u).QSimAnnPolBef(1)
                        % Save new solution because it shows better Q values
                        for lay = 2:Nlayer
                          
                            linkDQN(u).my_nnEv.wSv{lay} = linkDQN(u).my_nnEv.w{lay};  % Weights
                            linkDQN(u).my_nnEv.bSv{lay} = linkDQN(u).my_nnEv.b{lay};  % Biases
                        end    
                        for  xcnt = 1:NStates            
                            linkDQN(u).QSimAnnPolBef(xcnt) =  linkDQN(u).Qmat(xcnt,linkDQN(u).policyAsIndx(xcnt));
                        end    
                        linkDQN(u).QSimAnnPolAsInd = linkDQN(u).policyAsIndx;
                        linkDQN(u).QSimAnnPolicy = linkDQN(u).policy;
                    end

                    linkDQN(u).my_nnEv.a{1} = zeros( N_neuron_layer(1) , 1 );  
                    for lay = 2:Nlayer

                        linkDQN(u).my_nnEv.w{lay} = rand( N_neuron_layer(lay) , N_neuron_layer(lay-1) );  % Weights
                        linkDQN(u).my_nnEv.b{lay} = rand( N_neuron_layer(lay) , 1);  % Biases
                        linkDQN(u).my_nnEv.a{lay} = zeros( N_neuron_layer(lay) , 1 );  % Outputs (after activation function)
                    end        

                    linkDQN(u).QmatEntry = zeros(NStates,length(ppossible)); % Initialize visit counter in Q-tables
                    % Make forward pass to get vector of Q-values for each state. 
                    for  xcnt = 1:NStates
                        linkDQN(u).my_nnEv = make_nn_fwd_pass(xcnt,linkDQN(u).my_nnEv);
                        linkDQN(u).Qmat(xcnt,:) = linkDQN(u).my_nnEv.a{Nlayer};
                    end
                    linkDQN(u).QmatTrgt = linkDQN(u).Qmat;  % Update target Q-values
                    linkDQN(u).policy = ppossible(randi(length(ppossible),1,NStates)); 
                    [RowInd,PolIndxVal]=find(bsxfun(@eq,[ppossible;ppossible],linkDQN(u).policy'));  % AK 10/21/2019
                    linkDQN(u).policyAsIndx(RowInd) = PolIndxVal; % AK 10/21/2019
                    linkDQN(u).ActAvailbl = [ones(NStates,1) * ppossible , linkDQN(u).policy']; % Each row a different state

                    linkDQN(u).LearnRateVec = LearnRate * ones(length(ppossible),1);
                    xt=ones(1,NCRTxRx); % Initial state, with no CR transmission yet, the primary throughput can't have changed
                    xtsv = xt(1);
               end
            end

            % Load the best solution from "simulated annealing"
            if c == LenAllSimAnn 
                for u=1:NCRTxRx
                    for lay = 2:Nlayer
                       
                        linkDQN(u).my_nnEv.w{lay} = linkDQN(u).my_nnEv.wSv{lay};  % Weights
                        linkDQN(u).my_nnEv.b{lay} = linkDQN(u).my_nnEv.bSv{lay};  % Biases
                    end        
                    linkDQN(u).QmatEntry = zeros(NStates,length(ppossible)); % Initialize visit counter in Q-tables
                    % Make forward pass to get vector of Q-values for each state. 
                    for  xcnt = 1:NStates
                        linkDQN(u).my_nnEv = make_nn_fwd_pass(xcnt,linkDQN(u).my_nnEv);
                        linkDQN(u).Qmat(xcnt,:) = linkDQN(u).my_nnEv.a{Nlayer};
                    end
                    linkDQN(u).QmatTrgt = linkDQN(u).Qmat;  % Update target Q-values
                    linkDQN(u).policy = linkDQN(u).QSimAnnPolicy; 
   
                    linkDQN(u).policyAsIndx = linkDQN(u).QSimAnnPolAsInd;
                    linkDQN(u).ActAvailbl = [ones(NStates,1) * ppossible , linkDQN(u).policy']; % Each row a different state
  
                    linkDQN(u).LearnRateVec = LearnRate * ones(length(ppossible),1);
                    xt=ones(1,NCRTxRx); % Initial state, with no CR transmission yet, the primary throughput can't have changed
                    xtsv = xt(1);
                end
            end

            c = c + 1;
        end  % End of all exploration phases
            
        % Evaluate result deep Q-learning network
        for n = 1:NCRTxRx
            PTxResDQN(n) = linkDQN(n).policy(1); % Power at CRs
            ActionIndxCRFromDQN(mcrun,n) = linkDQN(n).policyAsIndx(1);
        end    
        PTxCR = pow2db(PTxResDQN) + 30; % Power at CRs converted to dBm 
        NoiseAtPrmy = PTxResDQN * CRtoPrimChannel' + db2pow(Pn-30); % Calculate noise at primary UEs, secondary networks appears as noise
        sinrUE = pow2db(  TotUERxTxP ./(TotUEIntf + NoiseAtPrmy) );
        % Get throughput and modulation index in primary network
        [ThrgptPrimDQN,IndxModTypeDQN] = GetThrghptandModType(sinrUE,SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,ModType_AMC_PedB_3retx_2x2_BLER1em1,CQI_AMC_PedB_3retx_2x2_BLER1em1);
        for n = 1:NCRTxRx
            IntfCh = CRChannel(n,:);
            IntfCh(n) = [];
            IntfPw = db2pow(PTxCR-30);
            IntfPw(n) = [];
            SINRDen(n) = IntfCh * IntfPw' + PrimtoCRChannel(n,:) * db2pow(PTx-30)' + db2pow(Pn-30); % Primary network appears as noise
        end
        SINRCRRxDQN = diag(CRChannel)' .* db2pow(PTxCR-30) ./ SINRDen;
        % Calculate throughput at CRs
        ThrgptCRStd = GetThrghptandModType(pow2db(SINRCRRxDQN),SNR_PedB_3retx_2x2,T_AMC_PedB_3retx_2x2_BLER1em1,ModType_AMC_PedB_3retx_2x2_BLER1em1,CQI_AMC_PedB_3retx_2x2_BLER1em1);

        RelThrChngDQN = abs(ThrgptPrimNoCR - ThrgptPrimDQN)./ ThrgptPrimNoCR;
        % record worst case of relative throughput case
        MAXRelThrChngDQN(mcrun) = max(RelThrChngDQN);   
        CRThrghFromDQN(mcrun,:) = ThrgptCRStd;
        RelThrChngNearestPNDQN(mcrun,:) = RelThrChngDQN(ClosestCRTx);

        
        % -------------------------------------------------------
        % -------------------------------------------------------
        % ++++++++++  End of Deep Q-learning Network  +++++++++++
        % -------------------------------------------------------
        % -------------------------------------------------------
            
        if FirstPass
            ThrgptPrimFirstPass(mcrun) = mean(ThrgptPrim);
            ThrgptCRFirstPass(mcrun) = mean(ThrgptCR);
            sinrUEFirspPass = [sinrUEFirspPass sinrUE]; 
            FirstPass = 0;
        end    

        CRTXPwrdBm(:,mcrun) = PTxCR';

        ThrgptRunWithCRTXonModT0(mcrun) = mean(ThrgptPrim);
        SINRPrimWithCRTXonModT0 = [SINRPrimWithCRTXonModT0 sinrUE];
        CRThrgptRunWithCRTXonModT0(mcrun) = mean(ThrgptCR);

        RelThrChngRun(mcrun,:) = RelThrChngPhase;
        disp('RelThrChngRun')
        RelThrChngRun(mcrun,end)

        for u = 1:NCRTxRx
            link(u).QmatTrackAppAvg(mcrun,:) = mean(link(u).QmatTrackApp);
            link(u).QmatTrackAppStd(mcrun,:) = std(link(u).QmatTrackApp);
        end
        
        % Wrap up
        ThrgptRunWithCR(mcrun) = mean(ThrgptPrim);
%        ThrgptRunWithCRFinal(mcrun) = mean(ThrgptPrimFinal);
        ThrgptRunAtCRFinal(mcrun) = mean(ThrgptCR);
        ThrgptRunAtCR(mcrun) = mean(ThrgptCRStd);
        SINRPrimWithCR = [SINRPrimWithCR sinrUE];
        [mcrun , IterCnt , ThrgptRunNoCR(mcrun) , ThrgptPrimFirstPass(mcrun) , ThrgptRunWithCRTXonModT0(mcrun), ThrgptRunWithCR(mcrun) , ThrgptCRFirstPass(mcrun) , CRThrgptRunWithCRTXonModT0(mcrun) ThrgptRunAtCR(mcrun)]
        disp('MAXRelThrChngDQN')
        MAXRelThrChngDQN(mcrun)

        CCatPolSv = [];
        for u = 1:NCRTxRx
            CCatPolSv = [CCatPolSv;linkDQN(u).PolicySv'];
            AvgPseudoSNRQval(mcrun,u) = mean(mean((10*log10(linkDQN(u).QmatSv(EpLen:2:end,:).^2./movstd(linkDQN(u).QmatSv(EpLen:2:end,:),EpLen).^2))));
            linkDQN(u).MovMeanRelLearnError(mcrun,:) = movmean(linkDQN(u).Relerror_epoch(mcrun,:),EpLen);
            meanRunPseudoPSNQvalPol(mcrun,u) = mean(10*log10(linkDQN(u).QmatSvPol(:,1).^2./movstd(linkDQN(u).QmatSvPol(:,1),EpLen).^2));
            SizeNextPolcySet(u+(mcrun-1)*NCRTxRx,:) = linkDQN(u).SizeNextPolcySetSv;
            sortQval = sort(linkDQN(u).QmatSv(1:2:end,:),2,'descend');
            SecndLargstQ(u+(mcrun-1)*NCRTxRx,:) = sortQval(:,2)';
            DiffQPolAndNextLgst(u+(mcrun-1)*NCRTxRx,:) = linkDQN(u).QmatSvPol(:,1)-SecndLargstQ(u+(mcrun-1)*NCRTxRx,:)'; % Difference in the Q-value between policy and next largest (interest to see if derivative goes to 0
        end
       ChPol = sum(CCatPolSv(:,2:end)-CCatPolSv(:,1:end-1));
        ChPol = sum(abs(CCatPolSv(:,2:end)-CCatPolSv(:,1:end-1)));  % AK: Fixed bug 11/9/19
        PolConvg(mcrun) = max(find(ChPol))+1;
        RelDifCRSumThrgh(mcrun) = (sum(BestThrgptEachCR(mcrun,:))-sum(CRThrghFromDQN(mcrun,:)))/sum(BestThrgptEachCR(mcrun,:)); %AK 10/21/2019


        disp('~~~~~~~~~~~~~~~~~~~~~~~~');
        disp('Rel. Throughput Change PN Q-table ');RelThrChngPhase(end)
        disp('CR Throughput from policy Q-table ');CRThrghFromPolcy(end,:)

        disp('Rel. Throughput Change PN DQN ');MAXRelThrChngDQN(end)
        disp('CR Throughput from policy DQN ');CRThrghFromDQN(end,:) 
        disp('Final CR policy (as index to ppossible) DQN ');ActionIndxCRFromDQN(end,:)
        disp('Rel. Throughput Change Nearest PN DQN ');RelThrChngNearestPNDQN

        figure(5);plot(1:length(linkDQN(1).QmatSv(1:2:end,:)),linkDQN(1).QmatSv(1:2:end,:),EpLen:EpLen:length(linkDQN(1).QmatSv(1:2:end,:))+1,linkDQN(1).TolrcBorder(1,:),'-o');legend;grid;title('Q-values, CR 1');
%         figure(6);semilogy(1:length(linkDQN(1).Relerror_epoch),linkDQN(1).Relerror_epoch,1:length(linkDQN(1).Relerror_epoch),movmean(linkDQN(1).Relerror_epoch,EpLen),'-r.',1:length(linkDQN(1).Relerror_epoch),movstd(linkDQN(1).Relerror_epoch,EpLen),'-k.');grid;legend
%         title('relative error, its moving mean and moving standard deviation, CR 1');
%         figure(7);plot(linkDQN(1).SaveRwrd');grid;legend;title('rewards, CR 1');
%         figure(8);plot(linkDQN(1).QmatSTD);grid;legend;title('Q values standard deviation, CR 1');
%         figure(9);plot(1:4800,((linkDQN(1).QmatSvPol(:,1)-SecndLargstQ(end-1,:)')),EpLen*5:EpLen*5:4800,gradient((linkDQN(1).QmatSvPol(1:EpLen*5:end,1)-SecndLargstQ(end-1,1:EpLen*5:end)')));grid
        figure(9);plot(10*log10(linkDQN(1).QmatSvPol(:,1).^2./movstd(linkDQN(1).QmatSvPol(:,1),EpLen).^2));legend;grid;title('Sort of Q-values SNR [dB], CR 1');
        figure(15);plot(1:length(linkDQN(2).QmatSv(1:2:end,:)),linkDQN(2).QmatSv(1:2:end,:),EpLen:EpLen:length(linkDQN(2).QmatSv(1:2:end,:))+1,linkDQN(2).TolrcBorder(1,:),'-o');legend;grid;title('Q-values, CR 2');
%         figure(16);semilogy(1:length(linkDQN(2).Relerror_epoch),linkDQN(2).Relerror_epoch,1:length(linkDQN(2).Relerror_epoch),movmean(linkDQN(2).Relerror_epoch,EpLen),'-r.',1:length(linkDQN(2).Relerror_epoch),movstd(linkDQN(2).Relerror_epoch,EpLen),'-k.');grid;legend
%         title('relative error, its moving mean and moving standard deviation, CR 2');
%         figure(17);plot(linkDQN(2).SaveRwrd');grid;legend;title('rewards, CR 2');
%         figure(18);plot(linkDQN(2).QmatSTD);grid;legend;title('Q values standard deviation, CR 2');
%        figure(19);plot(10*log10(linkDQN(2).QmatSv(1:2:end,:).^2./movstd(linkDQN(2).QmatSv(1:2:end,:),EpLen).^2));legend;grid;title('Sort of Q-values SNR [dB], CR 2');
%         figure(19);plot(1:4800,((linkDQN(2).QmatSvPol(:,1)-SecndLargstQ(end,:)')),EpLen*5:EpLen*5:4800,gradient((linkDQN(2).QmatSvPol(1:EpLen*5:end,1)-SecndLargstQ(end,1:EpLen*5:end)')));grid
        figure(19);plot(10*log10(linkDQN(2).QmatSvPol(:,1).^2./movstd(linkDQN(2).QmatSvPol(:,1),EpLen).^2));legend;grid;title('Sort of Q-values SNR [dB], CR 2');



%         figure(9);plot(DiffQPolAndNextLgst(106*2-1,:));grid
%         figure(19);plot(DiffQPolAndNextLgst(106*2,:));grid        
        
        for c = 1:length(ppossible)
            for u = 1:NCRTxRx
                meanRwrd(u,c) = mean( linkDQN(u).SaveRwrd(c , 1:linkDQN(u).SaveRwrdIndx(c)-1 ));
                StdRwrd(u,c) = std( linkDQN(u).SaveRwrd(c , 1:linkDQN(u).SaveRwrdIndx(c)-1 ));
            end
        end
        meanRwrd
        StdRwrd
        

       actionDQNindex_and_ToState = [actionDQNindexSv, xtsv'];  % Gives error, During exploration in DQN: actions taken (as index into ppossible) and resulting state
        disp('BestRelThrChng');
        
        BestRelThrChng
        disp('BestThrgptEachCR');
        BestThrgptEachCR
        disp('BestCRPwrIndxSv');
        BestCRPwrIndxSv
        disp('BestRelThrChngNearestPNExCase');
        BestRelThrChngNearestPNExCase
        
        [linkDQN(1).PolicySv';linkDQN(2).PolicySv']

        RelDifCRSumThrgh

        [ActionIndxCRFromDQN,BestCRPwrIndxSv,(ActionIndxCRFromDQN-BestCRPwrIndxSv),BestThrgptEachCR,sum(BestThrgptEachCR,2),CRThrghFromDQN,sum(CRThrghFromDQN,2),RelDifCRSumThrgh']

        if RelDifCRSumThrgh(mcrun) ~= 0
            RelDifCRSumThrgh(mcrun)
        end
        
    end


    toc

% end
for u = 1:NCRTxRx
    MCMean_MovMeanRelLearnError(u,:) = mean(linkDQN(u).MovMeanRelLearnError,1);
    meanPseudoPSNQvalPol(u) = mean(10*log10(linkDQN(u).QmatSvPol(:,1).^2./movstd(linkDQN(u).QmatSvPol(:,1),EpLen).^2));
end
figure();semilogy(MCMean_MovMeanRelLearnError');grid



disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');

% AvgPseudoSNRQval
meanRunPseudoPSNQvalPol
meanPseudoPSNQvalPol = mean(meanRunPseudoPSNQvalPol,1)
PolConvg
RelDifCRSumThrgh
% 
% disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
% diary off;
% save(WrkSPcSvFile);
% 
