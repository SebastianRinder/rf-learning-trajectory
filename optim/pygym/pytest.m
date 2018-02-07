clear classes;
P = py.sys.path;
% if count(P,'/home/yx58adif/openaigym') == 0
%     insert(P,int32(0),'/home/yx58adif/openaigym');
% end
% if count(P,'/home/yx58adif/openaigym/gym') == 0
%     insert(P,int32(0),'/home/yx58adif/openaigym/gym');
% end
% if count(P,'/home/yx58adif/imageio') == 0
%     insert(P,int32(0),'/home/yx58adif/imageio');
% end
% if count(P,'/home/yx58adif/mujoco-py') == 0
%     insert(P,int32(0),'/home/yx58adif/mujoco-py');
% end
% if count(P,'/home/yx58adif/Box2D-kengz-2.3.3') == 0
%     insert(P,int32(0),'/home/yx58adif/Box2D-kengz-2.3.3');
% end
% if count(P,'/home/yx58adif/pyglet') == 0
%     insert(P,int32(0),'/home/yx58adif/pyglet');
% end
% if count(P,'/home/yx58adif/six') == 0
%     insert(P,int32(0),'/home/yx58adif/six');
% end
% if count(P,'/home/yx58adif/numpy') == 0
%     insert(P,int32(0),'/home/yx58adif/numpy');
% end



disp(P);

mod = py.importlib.import_module('gymCartPole');
py.reload(mod);
py.gymCartPole.init_environment('CartPole-v0');
policy= [2, 3, 1, 0, 1, 1, 1, 0];
timesteps = uint16(1000);
actionlist = uint8([0,1]);
ret = py.gymCartPole.evaluate(policy,timesteps,actionlist);

states = double(py.array.array('d',py.numpy.nditer(ret{1,1})));
trajectory.state = reshape(states,4,[])';
trajectory.action = double(py.array.array('d',py.numpy.nditer(ret{1,2})))';
trajectory.prob = double(py.array.array('d',py.numpy.nditer(ret{1,3})))';
trajectory.cumReward = double(py.array.array('d',py.numpy.nditer(ret{1,4})))';
