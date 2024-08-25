module saveSettings

    export saveInRealUnits, saveData, saveDataFolderTag, addHeaders

    saveInRealUnits::Bool = true # if true, save vel+accel in m/s, mm/ms^2.  If false, save in normalized units (vel= v/(gam/k)), (force=1e-3*hbar*k*gam)
    saveData::Bool = true # if you want to save the data
    saveDataFolderTag::String = "SrFRedMOTNormalValues" # If you want to put anything additional in "savefoldername" to tag it, see variable folderString after lasers are declared.
    addHeaders::Bool = true

end