## Wrapper function for calculating different land surface parameter
lsp <- function(in_dem, 
                calculate_all = FALSE,
                hillshade = NULL,
                slope = NULL,
                twi = NULL,
                pisr = NULL, 
                tpi = NULL, 
                rel_ele = NULL,
                tpi_radius = 50,
                rel_ele_radius = 50,
                rm_twi_tmp_files = FALSE,
                use_sca = F,
                pisr_start_date = 20, 
                pisr_end_date = 23, 
                pisr_start_month = 3,
                pisr_end_month = 9,
                pisr_d_step = 6,
                pisr_t_step = 4,
                pisr_year = 2019,
                pisr_latitude = 69,
                log_file = "log_file.txt"
                ) 
  {
  # This function calculates desired lsp's from an input DEM (in_dem)
  # If an output filename is given the program calculates the specific LSP.
  # Check the parameters and usage for each variable separately from the 
  # documentation.
  # Default PISR parameters follow Kemppinen et al. 2018. 
  require(RSAGA)
  require(raster)
  
  try(if(file.exists(log_file)) stop("ERROR! Log file exists, give new name!"))
  
  sink(log_file)
  # Import digital elevation model if SAGA modules are being used (but not if only relative elevation is used)
  # Give default names if calculate_all is set to TRUE:
  if(calculate_all == TRUE){
    slope = "slope"
    twi = "twi"
    pisr = "pisr"
    hillshade = "hs"
    # Name TPI and relative elevation with given radius parameter (or by default radius if none is given)
    tpi = paste0("tpi_", tpi_radius)
    rel_ele = paste0("relative_elevation_", rel_ele_radius)
  }
  
  # Don't import DEM if only relative elevation is calculated (it uses raster-package)
  if( (!is.null(c(hillshade,slope,pisr,tpi,twi))) & !(file.exists("dem.sgrd")) ){
    rsaga.import.gdal(in_dem, "dem" , env=myenv)
  }
  
  ## Calculate hillshade
  if(!is.null(hillshade)){
    if(file.exists(paste0(hillshade, ".sgrd"))){
      print("Hillshade file already exists!")
    } else {
      rsaga.hillshade(in.dem = "dem.sgrd", out.grid = hillshade, env = myenv)
    }
  }
  
  ## Calculate slope
  # Slope reference:
  # Zevenbergen, L.W., Thorne, C.R., 1987. Quantitative-Analysis of Land Surface-Topography. Earth Surface Processes and Landforms 12, 47-56. https://doi.org/10.1002/esp.3290120107 .
  if(!is.null(slope)){
    if(file.exists(paste0(slope, ".sgrd"))){
      print("Slope file already exists!")
    } else {
    # Slope from the unfilled DEM
    rsaga.slope(in.dem = "dem.sgrd",
                out.slope = slope,
                method = "poly2zevenbergen",
                env = myenv
                )
    }
  }
  
  ## Calculate Topographic Position Index
  if(!is.null(tpi)){
      # if no radius is given default value of 50 m is used (prints warning)
      if(tpi_radius == 50){
        # print warning
        print("No radius was given! Using default value of 50 meters.")
      }
      if(file.exists(paste0(tpi, ".sgrd"))){
        print("TPI file already exists!")
      }else{
      # Radius is in map units.
      rsaga.get.usage("ta_morphometry", 18, env = myenv)
      rsaga.geoprocessor("ta_morphometry", 18, 
                     param = list(
                       DEM = "dem.sgrd",
                       TPI = tpi,
                       RADIUS_MIN = 0,
                       RADIUS_MAX = tpi_radius),
                     show.output.on.console = T, invisible = TRUE, intern = TRUE, env = myenv)
      }
  }
  
  ## Relative elevation
  # Create function (uses R and raster package)
  relative_elevation <- function(in_dem, radius, out_name) {
    require(raster)
    # Get resolution
    dem <- raster(in_dem)
    reso <- abs(dem@extent@xmin-dem@extent@xmax)/ncol(dem)
    
    # Calculate window size (2 x radius), units converted to pixels
    ws <- round(2*(radius/reso))
    
    # Force window size to odd number (3,5,7 ... )
    if((ws %% 2) == 0) { 
      ws <- ws+1 
    }
    
    # Calculate relative elevation
    minimum_elevation <- focal(dem, w=matrix(1, ws,ws), fun = min) # 1 = weight
    relative_elevation <- (dem-minimum_elevation)

    # Write out
    writeRaster(relative_elevation, filename = out_name, format="SAGA", overwrite=T)
  }
  
  # Run function:
  if(!is.null(rel_ele)){
    # if no radius is given default value of 50 m is used (prints warning)
    if(rel_ele_radius == 50){
      print("No radius was given! Using default value of 50 meters.")
    }
    if(file.exists(paste0(rel_ele, ".sgrd"))){
      print("Relative elevation file already exists!")
    } else {
      relative_elevation(in_dem, radius = rel_ele_radius, out_name = rel_ele)
    }
  }
  
  ## Calculate Topographic Wetness Index
  if(!is.null(twi)){
    if(file.exists(paste0(twi, ".sgrd"))){
      print("TWI file already exists!")
    } else {
    ## Step 1. Preprocess DEM
    #?rsaga.fill.sinks
    rsaga.fill.sinks(in.dem = "dem.sgrd",
                     out.dem = "dem_filled.sgrd",
                     method = "xxl.wang.liu.2006",
                     minslope = 0.01,
                     env = myenv
    )
    
    ## Step 2. Calculate slope (Zevenbergen & Thorne, 1987)
    # output in radians
    rsaga.slope(in.dem = "dem_filled.sgrd",
                out.slope = "slope_filled_radians.sgrd",
                method = "poly2zevenbergen",
                env = myenv
    )
    
    # Step 3. Calculate total catcment area, TCA  
    # with SAGA v. 2.1.3 or newer use rsaga.topdown.processing
    # convergence Follows Freeman 1991 (Quinn et al. would be 1.0)
    # linear.threshold # Disabled, see Sorensen et al. 2006 for discussion, (also Quinn et al., 1995 for background)
    rsaga.geoprocessor("ta_hydrology", 0,
                       param = list(
                         ELEVATION="dem_filled.sgrd",
                         FLOW="tca.sgrd",
                         METHOD=4,
                         LINEAR_DO=0,
                         CONVERGENCE=1.1,
                         NO_NEGATIVES=1
                       ), 
                       show.output.on.console = T, invisible = TRUE, intern = TRUE, 
                       env = myenv)
    
    if(use_sca == TRUE){
      # Step 4. Calculate specific catchment area (SCA), and flow width (optional)
      # TCA is input from step 3. 
      # for flow width calculation see Geomorphometry 2009, p. 181 (Gruber & Peckham, 2009)
      # method 1 = Quinn et al. 1991 The Prediction of Hillslope Flow Paths for Distributed Hydrological Modeling Using Digital Terrain Models. Hydrological Processes 5, 59-79. https://doi.org/10.1002/hyp.3360050106 
      rsaga.geoprocessor("ta_hydrology", 19,
                         param = list(
                           DEM="dem_filled.sgrd",
                           WIDTH="flow_width.sgrd",
                           TCA="tca.sgrd", 
                           SCA="sca.sgrd", 
                           METHOD=1
                         ), 
                         show.output.on.console = T, invisible = TRUE, intern = TRUE, 
                         env = myenv)
      
      # Step 5. Calculate Topographic Wetness Index
      # Riihimäki et al. paper used TCA/L ( = pSCA)
      # In practice there was little difference betweem pSCA and SCA in predictive performance (tested with maisemaikkuna-data)
      rsaga.get.usage("ta_hydrology",20, env=myenv)
      rsaga.geoprocessor("ta_hydrology", 20, 
                         param = list(
                           SLOPE="slope_filled_radians.sgrd",
                           AREA="sca.sgrd", 
                           TWI= twi, 
                           CONV=0,  
                           METHOD=0
                           ),
                         show.output.on.console = T, invisible = TRUE, intern = TRUE, env = myenv)
      } else {
      # Step 4. Calculate Topographic Wetness Index
      # Using TCA/L 
      # following Riihimäki et al. paper 
      # rsaga.get.usage("ta_hydrology",20, env=myenv)
      rsaga.geoprocessor("ta_hydrology", 20, 
                         param = list(
                           SLOPE="slope_filled_radians.sgrd",
                           AREA="tca.sgrd", 
                           TWI= twi, 
                           CONV=1,  
                           METHOD=0
                         ),
                         show.output.on.console = T, invisible = TRUE, intern = TRUE, env = myenv)
      }
    }
    if(rm_twi_tmp_files == T){
      # Remove temporary files assosiated with TWI processing
      file.remove(c("dem_filled.sgrd","dem_filled.sdat", "dem_filled.mgrd",
                    "slope_filled_radians.sgrd","slope_filled_radians.sdat", "slope_filled_radians.mgrd",
                    "flow_width.sgrd","flow_width.sdat", "flow_width.mgrd",
                    "sca.sgrd","sca.sdat", "sca.mgrd",
                    "tca.sgrd", "tca.sdat", "tca.mgrd")
        )
      }
    } # end of twi
  
    ## Calculate potential incoming solar radiation
    # Reference 
    if(!is.null(pisr)){
      if(file.exists(paste0(pisr, ".sgrd"))){
        print("PISR file already exists")
      }else{
      # See Bohner and Selige, 2009
      # Calculate sky-view-factor (follows settings in Aalto et al. 2017)
      rsaga.geoprocessor("ta_lighting", 3, 
                         param = list(
                          DEM = "dem.sgrd",
                          SVF = "svf.sgrd",
                          NDIRS = 16,
                          RADIUS = 10000
                         ),
                         show.output.on.console = T, invisible = TRUE, intern = TRUE, env = myenv)
        
      # Notify user about the used parameters  
      print(paste0("Calculating Potential Incoming Solar Radiation from Months ", pisr_start_month, " to ", pisr_end_month, ", with day step of ", pisr_d_step, " and time step of ",pisr_t_step, " hours"))
      print(paste0("Start and end dates are for the first and last month are ", pisr_start_date," ",pisr_end_date, ", respectively"))  
      
      # Calculate total potential incoming solar radiation
      # Note rsaga.pisr2 has a date bug in v. 1.3, so the geoprocessor instead.
      rsaga.geoprocessor("ta_lighting", 2, 
                         param = list(
                           GRD_DEM = "dem.sgrd",
                           GRD_SVF = "svf.sgrd",
                           UNITS = 0,
                           GRD_TOTAL = pisr,
                           LATITUDE = pisr_latitude,
                           DAY = paste0(pisr_start_month,"/",pisr_start_date,"/", pisr_year),
                           DAY_STOP = paste0(pisr_end_month,"/",pisr_end_date,"/", pisr_year),
                           DAYS_STEP= pisr_d_step,
                           PERIOD=2,
                           HOUR_RANGE_MIN=0,
                           HOUR_RANGE_MAX=24,
                           HOUR_STEP = pisr_t_step
                         ),
                         show.output.on.console = T, invisible = TRUE, intern = TRUE, env = myenv)
      
      }
    } # end of pisr
  sink()
} #end of lsp-function 
  