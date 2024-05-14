function [molecule] = def_molecule(moleculeName)

if moleculeName == "CaF"
    molecule.kA = 2*pi/606e-9;
    molecule.kB = 2*pi/531e-9;
    molecule.gam = 2*pi*8.3e6;
    molecule.mass = 59*1.67e-27;


elseif moleculeName == "CaOH"
    molecule.kA = 2*pi/626.4e-9;
    molecule.kB = 2*pi/555.2e-9;
    molecule.gam = 2*pi*6.4e6;
    molecule.mass = 57*1.67e-27;

elseif moleculeName == "MgF"
    molecule.kA = 2*pi/359.3e-9;
    molecule.kB = 2*pi/269e-9;
    molecule.gam = 2*pi*22e6;
    molecule.mass = 43*1.67e-27;

elseif moleculeName == "SrF"
    molecule.kA = 2*pi/663e-9;
    molecule.kB = 2*pi/579e-9;
    molecule.gam = 2*pi*6.63e6;
    molecule.mass = 107*1.67e-27;

elseif moleculeName == "SrOH"
    molecule.kA = 2*pi/688e-9;
    molecule.kB = 2*pi/611e-9;
    molecule.gam = 2*pi*7e6;
    molecule.mass = 105*1.67e-27;

elseif moleculeName == "Ag"
    molecule.kA = 2*pi/333e-9;
    molecule.kB = 2*pi/333e-9;
    molecule.gam = 2*pi*22e6;
    molecule.mass = 107*1.67e-27;

elseif moleculeName == "Au"
    molecule.kA = 2*pi/268e-9;
    molecule.kB = 2*pi/268e-9;
    molecule.gam = 2*pi*26.2e6;
    molecule.mass = 197*1.67e-27;



end


end