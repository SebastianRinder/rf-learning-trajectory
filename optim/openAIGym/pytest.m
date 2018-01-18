clear classes;
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
