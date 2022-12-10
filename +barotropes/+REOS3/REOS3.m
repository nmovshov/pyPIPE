classdef REOS3 < barotropes.Tabular
    %REOS3 Tabular adiabats from the Rostock group.

    %% Properties
    properties
        T_1bar
        X
        Y
    end
    
    %% The constructor
    % This function reads in data from file and assigns to P_vals and rho_vals.
    % Isentrope files from Nadine contain two header line followed by a 5-column
    % array so using dlmread is easiest. The columns are:
    % rho[g/cc]    T[K]    P[GPa]    u[kJ/g]    s[kJ/g/K]
    methods
        function obj = REOS3(filename)
            if nargin == 0, return, end % matlab likes us to allow empty calls
            
            % Read in raw data
            raw = dlmread(filename,'',2,0);
            
            % Post processing
            rho = raw(:,1)*1000; % 1st col is rho in g/cm^3
            P = raw(:,3)*1e9; % 3rd col is P in GPa
            
            % Assign to base class properties
            obj.P_vals = P;
            obj.rho_vals = rho;
            
            % (Optional) select interpolation options
            obj.interpolation_method = 'linear'; % doc interp1 for options
            obj.extrapolation_method = 'extrap'; % doc interp1 for options
            
            % That's it but you can assign some meta data if desired
            obj.meta.raw_table_in = filename;
        end
    end
end
