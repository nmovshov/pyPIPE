classdef SCVH < barotropes.Tabular
    %SCVH Tabular adiabats from the SCVH equations of state.

    %% Properties
    properties
        T_1bar
        X
        Y
    end
    
    %% The constructor
    % This function reads in data from file and assigns to P_vals and rho_vals.
    % Adiabat files from Chris contain one header line followed by a 3-column
    % array so using dlmread is very easy. However, the first non-header row
    % sometimes reports zero density at the maximum pressure, which will mess up
    % interpolation. So we skip that row. The last row is also apparently an
    % extrapolation to 1 bar from the original eos and sometimes fails. So we skip
    % the last line as well.
    % The columns are:
    % rho[g/cc]    P[dynes/cm^2]    T[K]
    methods
        function obj = SCVH(filename)
            if nargin == 0, return, end % matlab likes us to allow empty calls
            
            % Read in raw data
            raw = dlmread(filename,'',1,0);
            
            % Post processing
            rho = raw(2:end-1,1)*1000; % 1st col is rho in g/cm^3
            P = raw(2:end-1,2)/10;     % 2nd col is P in dynes/cm^2
            
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
