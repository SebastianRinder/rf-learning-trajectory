clear classes;
mod = py.importlib.import_module('gymCartPole');
py.reload(mod);
py.gymCartPole.init_environment('CartPole-v0');
py.gymCartPole.evaluate();