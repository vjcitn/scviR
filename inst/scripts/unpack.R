# this function retrieves key components from
# trained totalVI model for serialization; the totalVI
# references are invalid when session ends
unpack = function(mod, txbatch=c("PBMC10k", "PBMC5k")) {
 elbo_train = unname(unlist(mod$history_["elbo_train"]))
 elbo_val = unname(unlist(mod$history_["elbo_validation"]))
#
 latrep = mod$get_latent_representation()
 denoised = mod$get_normalized_expression(
    n_samples=25L, return_mean=TRUE, transform_batch=txbatch)
 names(denoised) = c("rna_denoised", "protein_denoised")
 protein_foreground = 100*model$get_protein_foreground_probability(
    n_samples=25L, return_mean=TRUE, transform_batch=txbatch)

 list(elbo_train=elbo_train, elbo_val=elbo_val, latrep=latrep,
     denoised=denoised, protein_foreground=protein_foreground)
}
