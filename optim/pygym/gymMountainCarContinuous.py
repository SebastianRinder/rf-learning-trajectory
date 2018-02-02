import numpy as np
import gym

env = []


def init_environment(str_environment):
    global env
    env = gym.make(str_environment)
    return env.action_space, env.observation_space


def evaluate(policy, timesteps, error_deviation):
    global env
    # actionlist = np.asarray(actionlist)
    policy = np.asarray(policy)
    state = reset()

    states = np.zeros((timesteps, np.shape(state)[0]))
    actions = np.zeros((timesteps, 1))
    probs = np.zeros((timesteps, 1))
    cum_rewards = np.zeros((timesteps, 1))

    cum_reward = 0

    for t in range(timesteps):
        #env.render()
        action, prob = action_selection_continuous(policy, state, error_deviation)
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


def action_selection_continuous(policy, state, error_deviation):
    feature = np.append(state, [1, state[0] ** 2, state[1] ** 2, state[0] * state[1], state[0] ** 2 * state[1], state[1] ** 2 * state[0], state[0] ** 3, state[1] ** 3])

    # error_deviation = 1e-1
    mu = np.dot(feature, policy)
    np.clip(mu, -1, 1)

    noise = np.random.normal(mu, error_deviation, 1) * error_deviation
    anext = mu + noise
    prob = -np.power(anext - mu, 2)  # / (2. * np.power(error_deviation, 2))
    return anext, prob


def reset():
    state = env.reset()
    return state

# init_environment('MountainCarContinuous-v0')
# evaluate(np.random.uniform(-1, 1, 10), 1000, np.array([0, 1, 2]))
