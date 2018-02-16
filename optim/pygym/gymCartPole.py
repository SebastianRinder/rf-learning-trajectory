import numpy as np
import gym

env = []


def init_environment(str_environment):
    global env
    env = gym.make(str_environment)
    return env.action_space, env.observation_space


def evaluate(policy, timesteps, actionlist):
    global env
    actionlist = np.asarray(actionlist)
    policy = np.asarray(policy)
    state = reset()

    states = np.zeros((timesteps, np.shape(state)[0]))
    actions = np.zeros((timesteps, 1))
    probs = np.zeros((timesteps, 1))
    cum_rewards = np.zeros((timesteps, 1))

    cum_reward = 0

    if isinstance(env.action_space, gym.spaces.discrete.Discrete):
        for t in range(timesteps):
            #env.render()
            action, prob = action_selection_discrete(policy, state, actionlist)
            state_next, reward, done, info = env.step(action)

            cum_reward = cum_reward + reward
            states[t, :] = state
            actions[t] = action
            probs[t] = prob
            cum_rewards[t] = cum_reward

            state = state_next
            if done:
                break


    return states[0:t, :], actions[0:t], probs[0:t], cum_rewards[0:t]


# def action_selection_continuous(policy, state, actionlist):
#     feature = np.concatenate((state, np.array([1])), axis=0)
#     policy = policy.reshape(actionlist.size, feature.size).T
#
#     error_deviation = 1e-4
#     mu = np.dot(state, policy) #np.transpose(policy)
#     mu[mu > 1] = 1
#     mu[mu < -1] = -1
#
#     noise = np.random.normal(mu, error_deviation, 1) * error_deviation
#     action_next = mu + noise
#     quad = -np.power(action_next - mu, 2) / (2. * error_deviation)


def action_selection_discrete(policy, state, actionlist):
    feature = np.concatenate((state, np.array([1])), axis=0)

    policy = policy.reshape(actionlist.size, feature.size).T
    p = np.exp(np.dot(feature, policy))
    p = np.divide(p, p.sum())

    tmp = np.random.rand() > np.cumsum(p)
    idx = tmp.sum()

    anext = actionlist[idx]
    prob = p[idx]
    return anext, prob


def reset():
    state = env.reset()
    return state

#init_environment('CartPole-v0')
#init_environment('MountainCar-v0')
#evaluate(np.array([2, 3, 1, 0, 1, -1,-2, 1, 1, 0]), 1000, np.array([0, 1]))
