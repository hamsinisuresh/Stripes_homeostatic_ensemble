# Stripes_homeostatic_ensemble
_Matlab scripts to simulate how entropic forces drive contact guidance in myofibroblasts_

cell_spread_shapes_metro1twoD_stripe() is the main function used to calculate distributions of cell morphometrics and free-energies for myofibroblasts on a stripe of given width 'swidth' (in microns).
Sample command to run this function: cell_spread_shapes_metro1twoD_stripe('width160', 0.2, 0.238, 1, 1e6, 160, 1.5, 1)


cell_spread_shapes_metro1twoD_stripe() calls on get_cell_energy2D_th_twoD_stripe() for expensive free-energy calculations. Mex get_cell_energy2D_th_twoD_stripe() to run this function faster. Paste lines 13-20 from get_cell_energy2D_th_twoD_stripe() on Matlab command prompt to generate mexed version of this function.

If mexing is not possible (for instance, when codegen is not installed/available), remove '_ mex' from lines 183 and 217 in cell_spread_shapes_metro1twoD_stripe() to run get_cell_energy2D_th_twoD_stripe() normally.


cytoskeletal_order_parameter() is used to assemble the overall distribution of stress-fibre orientations in cells on a given stripe from individual MCMC runs.
