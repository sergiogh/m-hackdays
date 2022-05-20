import json

input_file_name = "input.json"

with open(input_file_name) as f:
  dic = json.load(f)

if "extra_arguments" in dic:
    extra_arguments = dic['extra_arguments']
else:
    extra_arguments = {}


def build_QUBO(input_data, extra_arguments):
    # This is the core of your algorithm,
    # here you must return the corresponding answer.
    exp = input_data["exp"]
    cov = input_data["cov"]
    cost = input_data["cost"]
    budget = input_data["budget"]
    if ("theta0" in extra_arguments):
        theta0 = extra_arguments["theta0"]
    else:
        theta0 = 0.2
    if ("theta1" in extra_arguments):
        theta1 = extra_arguments["theta1"]
    else:
        theta1 = 0.4
    if ("theta2" in extra_arguments):
        theta2 = extra_arguments["theta2"]
    else:
        theta2 = 0.6

    n = len(exp)
    QUBO_def = {}

    for i in range(n):
        QUBO_def[(i, i)] = - exp[i] * theta0

    for i in range(n):
        QUBO_def[(i, i)] += cov[i][i] * theta1
        for j in range(i, n):
            QUBO_def[(i, j)] = 2 * cov[i][j] * theta1

    for i in range(n):
        QUBO_def[(i, i)] += (cost[i] ** 2 - 2 * budget * cost[i]) * theta2
        for j in range(i, n):
            QUBO_def[(i, j)] += 2 * cost[i] * cost[j] * theta2


    return QUBO_def

def run_sim(QUBO_def, provider="simmulated"):

    if provider == "hybrid":
        from dwave.system import LeapHybridSampler
        sampler = LeapHybridSampler()
    else:
        from neal import SimulatedAnnealingSampler
        sampler = SimulatedAnnealingSampler()

    sampleset = sampler.sample_qubo(QUBO_def)

    solution = []
    ans = sampleset.first.sample
    for k in ans:
        if ans[k] == 1:
            solution.append(k)

    return solution


QUBO = build_QUBO(dic['data'], extra_arguments)

result = run_sim(QUBO, "hybrid")
    

print(result)