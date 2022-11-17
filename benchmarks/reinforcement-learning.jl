using AscentOptimization
using BenchmarkTools
using Unitful: m, s 
model = Model(0.3/s,[0,2π])
explorer = EpsilonExplorer(0.1,0.0)
actions = AscentOptimization.action_space(model)
action=rand(actions)
vals = zeros(length(actions))
learner = AscentOptimization.QModel{AscentOptimization.RewardType}(actions)

@info "training model" @benchmark train_model!(model, standard_parameters, explorer)

@info "choice" @benchmark AscentOptimization.choice($explorer, $actions, $vals)

@info "randargmax" @benchmark AscentOptimization.randargmax($(rand(1:5, 10)))

@info "update!" @benchmark AscentOptimization.update!($learner, $action, 123m/s, 0.014)

@info "reward" @benchmark reward(model, standard_parameters)

@info "apply" @benchmark AscentOptimization.apply($model, $action)

