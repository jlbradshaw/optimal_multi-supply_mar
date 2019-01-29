# optimal_multi-supply_mar

Copyright Jonathan L. Bradshaw 2018

This is a repository for model scripts that describe and optimize designs for producing and delivering full advanced treated recycled water to existing stormwater spreading basins. runGeneralScenarios_paper2.m is the main script.

As a case study, the code uses the Hansen Spreading Grounds and three water recycling facilities in Los Angeles, CA.

The code is available here for your reference; unfortunately, I cannot provide support for its use.

This code was developed using MATLAB R2016b and relies on the Optimization Toolbox and Parallel Computing Toolbox.

Below is a summary of the purpose of each file:
•	computeLccFactor.m – Computes life cycle cost factors.
•	computeTotals.m – Aggregates data to facilitate analysis.
•	elevation_profile.m – Object with information about a pipeline route.
•	funit_lcc_daily_nonlinearCost.m – Function used for optimizing individual scenarios.
•	funit_lcc_daily_nonlinearCost_twoWRFparallel.m – Function used for optimizing parallel scenarios.
•	funit_lcc_daily_nonlinearCost_twoWRFserial.m – Function used for optimizing serial scenarios.
•	generalScenarios_paper2_par.m – Performs the system optimization process.
•	getConstraints_papertwo_singleWRF.m – Retrieves optimization constraints for individual systems.
•	getConstraints_papertwo_twoWRFparallel.m - Retrieves optimization constraints for parallel systems.
•	getConstraints_papertwo_twoWRFseries.m - Retrieves optimization constraints for serial systems.
•	getLccFactor.m – Retrieves life cycle cost factors.
•	hansen_sw_availability_afd_withDates_20171130.mat – Input data file with availability of stormwater in acre-feet per day.
•	hansen_sw_availability_cfs_withDates.mat – Input data file with availability of stormwater in cubic feet per second.
•	hansen_sw_availability_mgd_20171012.mat – Input data file with availability of stormwater in million gallons per day.
•	hessianfcn_daily.m – Function to compute Hessian for individual systems.
•	hessianfcn_dailytwoWRFparallel.m – Function to compute Hessian for individual systems.
•	hessianfcn_dailytwoWRFserial.m – Function to compute Hessian for individual systems.
•	import_piperoutes.m – Imports pipeline routes.
•	opt_paper2.m – Performs optimization routine for individual system.
•	opt_paper2_twoWRFparallel.m – Performs optimization routine for parallel systems.
•	opt_paper2_twoWRFseries.m – Performs optimization routine for serial systems.
•	pipe_op.m – Object with information about a pipeline.
•	pipeline_routine.m – Determines which pipe size options are feasible and compute the pumping requirements.
•	pipelineroutes.csv – Input data file containing pipeline routes to consider.
•	pipetemplate.csv – Input data file containing costs of piping.
•	prv.m – Object with information about a pressure release valve.
•	pumpdesign.m – Optimizes pipeline design.
•	pumpstation.m – Object with information about a pumpstation.
•	pumpstation_collection.m – Object with information about a set of pumpstations.
•	runGeneralScenarios_paper2.m – Main file.
•	units_check.m – Performs unit conversions.
