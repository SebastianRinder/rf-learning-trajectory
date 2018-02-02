import numpy as np
import gym

env = []


def init_environment(str_environment):
    global env
    env = gym.make(str_environment)


def evaluate(policy, timesteps, actionlist):
    global env
    actionlist = np.asarray(actionlist)
    policy = np.asarray(policy)
    state = reset()

    states = np.zeros((timesteps, np.shape(state)[0]))
    actions = np.zeros((timesteps, 4))
    probs = np.zeros((timesteps, 4))
    cum_rewards = np.zeros((timesteps, 1))

    cum_reward = 0

    for t in range(timesteps):
        env.render()
        action, prob = action_selection_continuous(policy, state, actionlist)
        state_next, reward, done, info = env.step(action)

        cum_reward = cum_reward + reward
        states[t, :] = state
        actions[t, :] = action
        probs[t, :] = prob
        cum_rewards[t] = cum_reward

        state = state_next
        if done:
            break

    return states[0:t, :], actions[0:t, :], probs[0:t], cum_rewards[0:t]


def action_selection_continuous(policy, state, actionlist):
    error_deviation = 1e-1
    #feature = np.concatenate((state, np.array([1])), axis=0)
    feature = np.append(state, 1)
    #feature = np.tile(feature, (4, 1))

    policy = policy.reshape(actionlist.size, feature.size).T

    mu = np.dot(feature, policy) #np.transpose(policy)
    mu[mu > 1] = 1
    mu[mu < -1] = -1

    noise = np.random.normal(0, 1, 4) * error_deviation
    anext = mu + noise
    prob = -np.power(anext - mu, 2) / (2. * np.power(error_deviation, 2))
    prob = np.exp(prob)/np.sqrt(2*np.pi*np.power(error_deviation, 2))
    return anext, prob


# def action_selection_discrete(policy, state, actionlist):
#     feature = np.concatenate((state, np.array([1])), axis=0)
#
#     policy = policy.reshape(actionlist.size, feature.size).T
#     p = np.exp(np.dot(feature, policy))
#     p = np.divide(p, p.sum())
#
#     tmp = np.random.rand() > np.cumsum(p)
#     idx = tmp.sum()
#
#     anext = actionlist[idx]
#     prob = p[idx]
#     return anext, prob


def reset():
    state = env.reset()
    return state

init_environment('BipedalWalker-v2')
evaluate(np.random.uniform(-1, 1, 100), 1000, np.array([0, 0, 0, 0]))
