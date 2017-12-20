import numpy as np
import gym

env = []


def init_environment(str_environment):
    global env
    env = gym.make(str_environment)


def evaluate(policy, timesteps, errorVariance):
    global env

    for t in range(timesteps):
        env.render()
        action = env.action_space.sample()
        state_next, reward, done, info = env.step(action)
        if done:
            break
        print(state_next)
        print(reward)
        print(done)
        print(info)


def action_selection_continuous(policy, state, action, errorvariance):
    mu = state
    np.dot(state, policy) #np.transpose(policy)
    mu[mu > 1] = 1
    mu[mu < -1] = -1

    noise = randn * errorvariance
    action_next = mu + noise
    quad = -np.power(action_next - mu, 2) / (2. * errorvariance)

def reset():
    state = env.reset()
    return state

init_environment('CartPole-v0')
state = reset()
evaluate()
