% +BAROTROPES
% A collection of pressure/density relations suitable for use in PIPE
% applications. Any subclass of +barotropes.Barotrope can be assigned to the eos
% field of TOFPlanet or CMSPlanet objects, as well as the fgeos and bgeos fields.
%
% Simple toy barotropes:
%   Barotrope        - Interface base class for all barotropes.
%   ConstDensity     - A constant density barotrope.
%   ConstPressure    - A constant pressure barotrope.
%   eos_table        - Example of a tabular barotrope that constructs from file.
%   IdealGasIsotherm - A toy barotrope - ideal gas at some constant T.
%   Polytrope        - A barotrope of the form P = K*rho^(1 + 1/n).
%   Tabular          - A base class for a generic table barotrope.
%
% +REOS3 subpackage with Rostock equation of state adiabats:
%   REOS3            - Subclass of barotropes.Tabular for REOS format tables.
%   REOS3.HHE        - Constructors for H/He REOS objects from table adiabats.
%
% +SCVH  subpackage with SCVH equation of state adiabats:
%   SCVH                 - Subclass of barotropes.Tabular for SCVH format tables.
%   SCVH.HHE             - Constructors for H/He SCVH objects from table adiabats.
%   SCVH.HE              - Constructors for SCVH Helium density on SCVH H/He adiabats.
%   SCVH.ANEOSSERPENTINE - Constructors for ANEOS serpentine density on SCVH H/He adiabats.
%   SCVH.ANEOSWATER      - Constructors for ANEOS water density on SCVH H/He adiabats.
%   SCVH.ANEOSIRON       - Constructors for ANEOS iron density on SCVH H/He adiabats.
