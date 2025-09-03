#Aim: get metrics from the .obj files through the habtools R package
#output file is called "results.txt"

##--------------------------------------------------
## Set working directory
setwd("E:/PhD_Organized/I. Science/E. Morphology/1. Models_obj_corrected")

library(Morpho)
library(raster)
library(habtools)
library(rgl)  # Needed for reading .obj files

mesh <- file2mesh("C1_2_t0.obj")

#check mesh resolution
my_fd <- function(file_paths) {
fractal_dimension = numeric(0)
resolution = numeric(0)
for (file in file_paths) {
  mesh <- file2mesh(file)
  print(basename(file))
  resvec <- Rvcg::vcgMeshres(mesh)[[2]]
  res <- mean(resvec[2])
  resolution <- rbind(resolution,
  file_name = basename(file), res) }
  return(resolution)
  
}

new_resolution <- my_fd(file_paths)
write.table(new_resolution, "E:/PhD_Organized/I. Science/E. Morphology/2. Traits_extraction/resolutions.txt",
            sep="\t", row.names = F)

file_paths <- list.files("E:/PhD_Organized/I. Science/E. Morphology/1. Models_obj_corrected", pattern = "\\.obj$", full.names = TRUE)

# Function to compute the desired metrics for each mesh file
compute_mesh_metrics <- function(file_paths) {
  
  # Create an empty dataframe to store results
  results <- data.frame(
    Sample_ID = character(0),
    planar_area = numeric(0),
    convexity = numeric(0),
    sphericity = numeric(0),
    packing = numeric(0),
    csf = numeric(0),
    fractal_dimension = numeric(0),
    surface_area = numeric(0),
    height_range = numeric(0),
    top_heaviness = numeric(0),
    rugosity = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each file path
  for (file in file_paths) {
    # Read the mesh file
    mesh <- readOBJ(file)
    
    # Compute planar area
    area <- planar(mesh)
    
    # Compute convexity
    convex <- convexity(mesh)
    
    # Compute sphericity
    sph <- sphericity(mesh)
    
    # Compute packing
    packing_value <- packing(mesh)
    
    # Compute CSF
    csf_value <- csf(mesh)
    
    # Compute fractal dimension
    fractal_dim <- fd(mesh, method = "cubes", diagnose = TRUE)
    # Compute surface area
    surf_area <- surface_area(mesh)
    
    # Apply meshtoDEM() function to the mesh
    dem <- mesh_to_dem(mesh)
    
    # Compute height range
    height_range_value <- hr(dem)
    
    # Compute top_heaviness
    top_h_value <-sma(mesh)
    
    # Compute rugosity
    rugosity_value <- rg(dem)
    
    # Add the result to the dataframe
    results <- rbind(results, data.frame(
      file_name = basename(file),
      planar_area = area,
      convexity = convex,
      sphericity = sph,
      packing = packing_value,
      csf = csf_value,
      fractal_dimension = fractal_dim,
      surface_area = surf_area,
      height_range = height_range_value,
      top_heaviness = top_h_value,
      rugosity = rugosity_value
    ))
  }
  
  # Return the dataframe
  return(results)
}

# Example usage
# Assuming your .obj files are in the "mesh_files" directory
file_paths <- list.files("E:/PhD_Organized/I. Science/E. Morphology/1. Models_obj_corrected", pattern = "\\.obj$", full.names = TRUE)
result_df <- compute_mesh_metrics(file_paths)

# View the result
print(result_df)
write.table(result_df, "E:/PhD_Organized/I. Science/E. Morphology/2. Traits_extraction/results_ok.txt",
            sep="\t", row.names = F)
