# core/integration.R
# Seurat layer integration and post-integration clustering workflow.

library(Seurat)
library(futile.logger)

# Integrate layers using a specified integration method, then find neighbors,
# clusters, and run UMAP. Retries integration with different k.weight values
# if the default fails.
integrateAndClusterSeurat <- function(
  seuratObject,
  integrationMethod,
  orig.reduction = "pca",
  new.reduction = "integrated",
  dims = 1:30,
  resolution = 0.5
) {
  flog.info("🚀 Starting integration and clustering.")

  k.weight.try <- c(NA, 50, 30)  # NA = default

  tryCatch({

    # Step 1: Integrate Layers (with retry)
    integration_success <- FALSE

    for (k in k.weight.try) {
      flog.info("🔄 Attempting integration. k.weight = %s",
                ifelse(is.na(k), "default", k))

      result <- tryCatch({
        if (is.na(k)) {
          IntegrateLayers(
            object = seuratObject,
            method = integrationMethod,
            orig.reduction = orig.reduction,
            new.reduction = new.reduction,
            verbose = TRUE
          )
        } else {
          IntegrateLayers(
            object = seuratObject,
            method = integrationMethod,
            orig.reduction = orig.reduction,
            new.reduction = new.reduction,
            k.weight = k,
            verbose = TRUE
          )
        }
      }, error = function(e) {
        flog.warn("⚠️ Integration failed with k.weight = %s : %s",
                  ifelse(is.na(k), "default", k),
                  e$message)
        return(NULL)
      })

      if (!is.null(result)) {
        seuratObject <- result
        integration_success <- TRUE
        flog.info("✅ Integration successful with k.weight = %s",
                  ifelse(is.na(k), "default", k))
        break
      }
    }

    if (!integration_success) {
      stop("Integration failed for all k.weight attempts.")
    }

    # Step 2: Re-join layers
    flog.info("🔗 Re-joining RNA layers.")
    seuratObject[["RNA"]] <- JoinLayers(seuratObject[["RNA"]])

    # Step 3: Find Neighbors
    flog.info("🔎 Finding neighbors using dims: %s", paste(dims, collapse = ", "))
    seuratObject <- FindNeighbors(seuratObject, reduction = new.reduction, dims = dims)

    # Step 4: Find Clusters
    flog.info("🧩 Finding clusters. Resolution = %s", resolution)
    seuratObject <- FindClusters(seuratObject, resolution = resolution)

    # Step 5: Run UMAP
    flog.info("🗺 Running UMAP.")
    if (length(dims) < 2) {
      flog.warn("⚠️ dims < 2 detected. Setting n.components = 1")
      seuratObject <- RunUMAP(seuratObject, dims = dims, reduction = new.reduction, n.components = 1)
    } else {
      seuratObject <- RunUMAP(seuratObject, dims = dims, reduction = new.reduction)
    }

    flog.info("🎉 Integration and clustering completed successfully.")
    return(seuratObject)

  }, error = function(e) {
    flog.error("❌ Fatal error: %s", e$message)
    stop("The integrateAndClusterSeurat function failed after retry attempts.")
  })
}
