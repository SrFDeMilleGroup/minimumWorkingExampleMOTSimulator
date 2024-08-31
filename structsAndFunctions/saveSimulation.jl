module saveSimulation

    using DelimitedFiles: writedlm
    using Dates: Dates, now

    using ..structs: Molecule, Lasers, GeneralSettings


    export saveToCsv


    function saveToCsv(mol::Molecule, lasers::Lasers, general::GeneralSettings, 
                        displacements_list::Vector{Float64}, userSpeeds_list::Vector{Float64}, longSpeeds_list::Vector{Float64},
                        obeResults::Vector{NTuple{11, Float64}})

        forceProjOnVelAvg = [result[1] for result in obeResults]
        forceProjOnVelUnc = [result[2] for result in obeResults]
        forceProjOnPosAvg = [result[3] for result in obeResults]
        forceProjOnPosUnc = [result[4] for result in obeResults]
        forceProjOnLongAvg = [result[5] for result in obeResults]
        forceProjOnLongUnc = [result[6] for result in obeResults]

        pExcAvg = [result[7] for result in obeResults]
        pF1DownAvg = [result[8] for result in obeResults]
        pF0Avg = [result[9] for result in obeResults]
        pF1UpAvg = [result[10] for result in obeResults]
        pF2Avg = [result[11] for result in obeResults]

        # create folder to save data
        folderString = string("./savedData/", general.simulationType, "bFieldSetting", general.bFieldSetting, 
                               "Force", general.forceProfile, "Dir", general.velDirRelToR, 
                               "NumLasers", lasers.numLasers, "Date", Dates.format(now(),"yyyymmdd_HHMMSS"))
        mkpath(folderString)

        # save OBE results
        open(string(folderString, "/obeResults", ".csv"), "a") do io
            headers = ["displacement (mm)",
                        "userSpeed (Gamma/kA)", "userSpeed (m/s)",
                        "longSpeed (Gamma/kA)", "longSpeed (m/s)",
                        "forceProjOnVelAve (1e-3\\hbar*k*\\Gamma)", "forceProjOnVelStd (1e-3\\hbar*k*\\Gamma)", 
                        "accelProjOnVelAve (mm/ms^2)", "accelProjOnVelStd (mm/ms^2)",
                        "forceProjOnPosAve (1e-3\\hbar*k*\\Gamma)", "forceProjOnPosStd (1e-3\\hbar*k*\\Gamma)",
                        "accelProjOnPosAve (mm/ms^2)", "accelProjOnPosStd (mm/ms^2)",
                        "forceProjOnLongAve (1e-3\\hbar*k*\\Gamma)", "forceProjOnLongStd (1e-3\\hbar*k*\\Gamma)",
                        "accelProjOnLongAve (mm/ms^2)", "accelProjOnLongStd (mm/ms^2)",
                        "PF1Down", "PF0", "PF1Up", "PF2", "PExc"]
            writedlm(io, reshape(headers, 1, :), ',') # make it a row, instead of a column

            writedlm(io, 
                    hcat(displacements_list, 
                        userSpeeds_list, userSpeeds_list .* mol.velFactor, 
                        longSpeeds_list, longSpeeds_list .* mol.velFactor, 
                        forceProjOnVelAvg, forceProjOnVelUnc, 
                        forceProjOnVelAvg .* mol.accelFactor, forceProjOnVelUnc .* mol.accelFactor, 
                        forceProjOnPosAvg, forceProjOnPosUnc, 
                        forceProjOnPosAvg .* mol.accelFactor, forceProjOnPosUnc .* mol.accelFactor, 
                        forceProjOnLongAvg, forceProjOnLongUnc, 
                        forceProjOnLongAvg .* mol.accelFactor, forceProjOnLongUnc .* mol.accelFactor, 
                        pF1DownAvg, pF0Avg, pF1UpAvg, pF2Avg, pExcAvg
                        ), 
                    ','
                    )
        end

        # save OBE settings
        open(string(folderString, "/obeSettings.csv"), "a") do io
            lasers_fields = fieldnames(Lasers)
            lasers_valuies = [string(getfield(lasers, field)) for field in lasers_fields]
            writedlm(io, ["Laser settings:"], ',')
            writedlm(io, [lasers_fields, lasers_valuies], ',')
            writedlm(io, "\n", ',')

            general_fields = fieldnames(GeneralSettings)
            general_values = [string(getfield(general, field)) for field in general_fields]
            writedlm(io, ["General settings:"], ',')
            writedlm(io, [general_fields, general_values], ',')
            writedlm(io, "\n", ',')

            mol_fields = fieldnames(Molecule)
            mol_values = [string(getfield(mol, field)) for field in mol_fields]
            writedlm(io, ["Molecule settings:"], ',')
            writedlm(io, [mol_fields, mol_values], ',')
        end
    end
end
