"""
Turbulence in velocity fields: bulk flows, structure functions, forward models, summaries and
inference. Every function takes ``rbcodes.ifu.vmap.VelocityMap`` objects (data or mocks,
image or source plane).

Modules
-------
bulk   bulk-flow models (Constant, Plane, Quadratic, ArctanDisk, TwoLobe) and compare_bulk
vsf    StructureFunction: raw / debiased / ivw / cross S2, jackknife covariance, directional
       S2, power-law fits (optionally with a noise floor or PSF smoothing), r_geo, r_turn
"""
