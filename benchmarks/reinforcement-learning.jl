using OrbitalMechanics
using BenchmarkTools

model = Model(0.3,[0,2π])
explorer = EpsilonExplorer(0.1,0.0)
actions = OrbitalMechanics.action_space(model)
action=rand(actions)
vals = zeros(length(actions))
learner = OrbitalMechanics.QModel(actions)

@info "training model" @benchmark train_model!(model, standard_parameters, explorer)

@info "choice" @benchmark OrbitalMechanics.choice($explorer, $actions, $vals)

@info "randargmax" @benchmark OrbitalMechanics.randargmax($(rand(1:5, 10)))

@info "update!" @benchmark OrbitalMechanics.update!($learner, $action, 123, 0.014)

@info "reward" @benchmark reward(model, standard_parameters)

@info "apply" @benchmark OrbitalMechanics.apply($model, $action)

