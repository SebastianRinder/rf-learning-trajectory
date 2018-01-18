% 	goto_pos       = 0,
% 	set_cannon_trg = 1,
% 	set_cannon_pos = 2,
% 	get_ball_state = 3,
% 	get_IK_hitPos  = 4,
% 	run_sim        = 5,
% 	set_ball_state = 6,
% 	execute_gains  = 7

% send gains to SL simulator
clear variables;

load('/home/akrour/postdoc/code/policysearchtoolbox/+Experiments/data/TableTennis_ClosedForm_varContext/TableTennis_TimeDependentREPS/settings001/eval001/trial001/trial.mat');
trial.GainProvider.applyNoise = 0;
trial.sampler.numSamples = 150;
% trial.contextSampler.BallNoiseAtBounce = .25;
% trial.contextSampler.MinBallPos = -.35;
% trial.contextSampler.MaxBallPos = -.35;
global ballAtTableHeightAfterStrikeIfAnyElseTableHeightAfterBounce;
ballAtTableHeightAfterStrikeIfAnyElseTableHeightAfterBounce = [];
data = trial.dataManager.getDataObject(0);
trial.sampler.createSamples(data);
rewards = data.getDataEntry('SLreturns');
plot (1:trial.sampler.numSamples, rewards, '.r');
plot(ballAtTableHeightAfterStrikeIfAnyElseTableHeightAfterBounce(:, 1), ...
    ballAtTableHeightAfterStrikeIfAnyElseTableHeightAfterBounce(:, 2), '.g');
hold on;
line([-0.763 0.763], [-3.64 -3.64]);
line([-0.763 -0.763], [-3.64 -.9]);
line([0.763 0.763], [-3.64 -.9]);
line([-0.763 0.763], [-.9 -.9]);
plot(0, -3.55, 'bo');
meanBallPos = mean(ballAtTableHeightAfterStrikeIfAnyElseTableHeightAfterBounce);
plot(meanBallPos(1), meanBallPos(2), 'rx');
line([-0.8 0.8], [-2.27 -2.27]);
hold off;
