
# plot_puck plots continuous or discrete values on a slide-seq puck

# (easter egg: If using categotical values, can include one category as '00NAN'
# for background grey points. This allows the broader puck orientation to be
# understood while focusing on a few cell types)

# val = value to plot. If factor -> discrete. If numerical -> numerical
# puck = seurat slide-seq puck with coords at thisPuck@images$image@coordinates. See [1] at end for how to create. If don't have this, look at coords argument
# barcodes = barcodes to plot over (if not include then all beads)
# pals = palette to plot. Can either pass a list of colors (for categorical) or a pals:: function (e.g pals::tableau20 , pals::kovesi.rainbow)
# alpha = alpha value for all beads
# raster = uses a nice package ggrastr to rasterize the puck output (so not thousands of points) but leave the text crisper when saving as pdf. Requires https://github.com/VPetukhov/ggrastr
# coords = alternative to passing in puck, just pass in coords
# alpha_bg = if include 00NAN background beads, their opacity

# requires library(pals), library(purrr)
# optionally:ggrastr (https://github.com/VPetukhov/ggrastr)
plot_puck <- function(val, puck=NULL, barcodes = NULL,
                      pal = pals::tableau20, # pal = pals::kovesi.rainbow,
                      alpha=0.001, raster=F, coords=NULL, alpha_bg = 0.2){
  # Take puck if available, then coords, or errors
  if(is.null(puck) & is.null(coords)){
    stop("NEED EITHER PUCK OR COORDS")
  }else if(!is.null(puck) ){
    my_table = puck@images$image@coordinates
  }else{
    my_table = coords
  }
  
  if(!is.null(barcodes)){
    my_table = my_table[barcodes,]
  }
  
  my_table$val = val
  
  # categorical value
  if(class(val) == "factor"){
    # equivalent to my_table$val = as.character(my_table$val)
    my_table$val %<>% as.character()
    
    # The legend is "CATEGORY_1 (10 beads)"
    tabled_vals = table(my_table$val)
    for(tbled_vals_nm in names(tabled_vals)){
      my_table[my_table$val == tbled_vals_nm, "val"] = glue("{tbled_vals_nm} ({tabled_vals[tbled_vals_nm]} beads)")
    }
    
    # Uses the sorted names of the categories. Change here if want different order
    my_table$val = factor(my_table$val, levels = sort(unique(my_table$val)))
    
    
    # If function then assume is pals:: function so call on how many categories = factor levels
    if(class(pal) == "function"){
      my_pal = pal(nlevels(my_table$val))
    }else{
      # If a list  then make sure the number of colors = number levels
      if(length(pal)!=nlevels(my_table$val)){
        stop(glue("NEED {nlevels(my_table$val)} manual classes in manual pal"))
      }
      my_pal = pal
    }
    
    print(length(levels(my_table$val)))
    
    # set names of palette to be the categories so ggplot knows how to map
    names(my_pal) = levels(my_table$val)
    # order palete to be the same as the categorical variables so legend matches up
    my_pal = my_pal[order(names(my_pal))]
    
    # If have a special category names 00NAN then is just background for orientation
    # So force those to be a light grey=#AAAAAA
    my_pal[grepl("00NAN", names(my_pal))] = "#AAAAAA"
    
    # alpha variable for beads
    my_table[["alpha"]] = alpha
    
    my_table[grepl("00NAN", my_table$val), "alpha"]=alpha_bg
    
    # CRITICAL: set order of points plotting (by alphabetical category name). By
    # 00NAN starting with 0, it goes first so the other categories can overlay
    # on top of it. Otherwise the background grey points can occlude
    
    point_order = seq(length(my_pal))
    names(point_order) = names(my_pal)
    
    my_table$order = point_order[my_table$val]
    my_table %<>% arrange(., order)
    
    # Easteregg, if want to force some categories to be on top/last in the
    # legend, then prefix the name with ZZ/. This will make it top and then
    # remove the ZZ/ prefix
    levels(my_table$val) %<>% str_remove("ZZ/")
    names(my_pal) %<>% str_remove("ZZ/")
    
    
    plot_scale = list(
      ggplot2::scale_color_manual(values=my_pal,
                                  limits=names(my_pal),
                                  breaks=discard(names(my_pal), ~grepl("00NAN", .))
      ))
    # If want to force the legend to be 3 rows for consistency, add this to the
    # above list
    # guides(color = guide_legend(nrow = 3)))
    
  }else{
    # If the value isn't a factor, it's continuous so plot differently
    
    print("cont")
    
    # CRITICAL: arrange value so high values are on top
    my_table %<>% arrange(., val)
    # Hack: if using the default tableau20 categorical color scheme, translate
    # to a default viridis here
    if(identical(pal, pals::tableau20)){
      pal = "viridis"
      plot_scale = ggplot2::scale_color_continuous(type =pal)
    }else{
      print("custom")
      plot_scale = ggplot2::scale_color_gradientn(colours = pal)
    }
  }
  
  # If raster=T and ggrastr installed then can use rasterized plotting
  point_fn = (ggplot2::geom_point)
  if(raster==T){
    if(!"ggrastr" %in% rownames(installed.packages())){
      stop("Need ggrastr package. Install with install.packages('ggrastr') or set raster=F")
    }
    point_fn = (ggrastr::geom_point_rast)
  }
  
  print (class(my_table))
  print (class(val))
  
  ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_point( ggplot2::aes(color=val), size=8.0, shape=19, show.legend = FALSE) +
    #scale_alpha_continuous(range = c(0.1, 0.2))
    scale_colour_gradient(
      
      low = "#F7F7F788",
      high = "#000000FF",
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour"
    ) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_blank(), 
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(), axis.title.x = element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(), axis.title.y = element_blank())
  
    
  # point_fn)
  #   ggplot2::aes(size = 1, shape=19, color=val, alpha=alpha))
  #raster.dpi=300,
  #bins=100) +
  # plot_scale+
  # ggplot2::scale_shape_identity() +
  # ggplot2::theme_classic() +
  # ggplot2::scale_size_identity() +
  # # Critical: if changing opacities otherwise weirdly tries to scale opacity
  # scale_alpha_identity() +
  # coord_fixed() +
  # theme(legend.position="top")+
  # labs(x = NULL, y = NULL)+
  # scale_x_continuous(expand=c(0,0)) +
  # scale_y_continuous(expand=c(0,0)) +
  # guides(x = "none", y = "none")+
  # # guides(alpha = "none")+
  # labs(color="Sub-cluster")
  
  
#   ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) +
#     point_fn(
#       ggplot2::aes(size = 1, shape=19, color=val, alpha=alpha),
#       raster.dpi=300,
#       bins=100) +
#     plot_scale+
#     ggplot2::scale_shape_identity() +
#     ggplot2::theme_classic() +
#     ggplot2::scale_size_identity() +
#     # Critical: if changing opacities otherwise weirdly tries to scale opacity
#     scale_alpha_identity() +
#     coord_fixed() +
#     theme(legend.position="top")+
#     labs(x = NULL, y = NULL)+
#     scale_x_continuous(expand=c(0,0)) +
#     scale_y_continuous(expand=c(0,0)) +
#     guides(x = "none", y = "none")+
#     # guides(alpha = "none")+
#     labs(color="Sub-cluster")
}