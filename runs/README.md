# Descriptions of Directories

I know the file naming convention isn't great, so I'm going to explain it here.

**Resolution:** The naming convention is `cluster_cooling_X` for X particles per 400 pc. (E.g. `cluster_cooling_200` has a resolution of 2 pc; 400 pc / 200 particles).

**Box Size:** Check `{run}/inputs/{run}.params.base`; generally the box size is 600 pc, except for "`large`" runs which are 1200 pc.

**"`small_steps`":** These are subsets of the normal runs, with higher time resolution for their outputs for a portion of time (and then not run to completion). This is generally to see behavior soon after a specific SN.

**"`early`":** Same as `small_steps` but for the 2nd SN. If the name has no "`early`", then it's probably after the 6th SN.

**"`cooling`":** Denotes whether radiative cooling was enabled. If the name doesn't include "`cooling`" then it's assumed to be adiabatic.

**`single`|`double`:** Test runs with just 1 or 2 SNe. Not used in the paper, but left in as a starting place for an early diagnostic

If you're still confused, feel free to email me and I can try to help.

