# using Revise package in mainSimulationCode.jl to monitor and update the changes in this script, to reduce the need to restart the kernel when making changes
# https://timholy.github.io/Revise.jl/stable/config/#Configuring-the-revise-mode
__revise_mode__ = :eval 

module generalSettingsOne

export bGradReal, waistInMM, numTrialsPerValueSet, velDirRelToR, initDispDir

bGradReal::Float64 = 12.5 # in units Gauss/cm. If bFieldSetting == "Static", this becomes the static field in Gauss
waistInMM::Float64 = 7 # mm, only used if polType is 3D. Handles finite MOT beam waists
numTrialsPerValueSet::Int64 = 100 # number of trials per set of values (displacementsInMM, userSpeeds, longSpeeds)
velDirRelToR::String = "Same" # relative direction of molecule velocity, w.r.t. initial displacement r. Options are ["Same", "Orthogonal", "Opposite", "Random"]
initDispDir::String = "XY" # if "XY", will force initial r to go along (x+y)/sqrt(2). Simulates slowing/trapping of molecules moving along slowing axis in tandem with velDirToR = "Same". If "Z", force go along Z

@assert waistInMM > 0 "Invalid waistInMM value: $waistInMM. It must be greater than 0."
@assert numTrialsPerValueSet > 0 "Invalid numTrialsPerValueSet value: $numTrialsPerValueSet. It must be greater than 0."
@assert velDirRelToR in ["Same", "Orthogonal", "Opposite", "Random"] "Invalid velDirRelToR value: $velDirRelToR. It must be one of ['Same', 'Orthogonal', 'Opposite', 'Random']."
@assert initDispDir in ["XY", "Z"] "Invalid initDispDir value: $initDispDir. It must be either 'XY' or 'Z'. "

end