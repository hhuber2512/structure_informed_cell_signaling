using NBInclude

#for analysis needing pairs (augmented and baseline comparison)
augmented_run = "202"
baseline_run = "203"
prior_nonbinding_offset = [3,3]
prior_binding_offset = [1,1];
@nbinclude("600_analysis.ipynb")

#run = 205
#@nbinclude("300_process_posterior_samples.ipynb")
#@nbinclude("500_generate_predictions.ipynb")
