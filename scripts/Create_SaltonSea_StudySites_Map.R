# Load necessary packages
library(mapview)
library(sp)
library(sf)
library(webshot2)

# NOTE: you need lat and long coordinates for this script!
# Create example data of points

#load("data/Metagenomes/SSD_mgm_analysis_data.Rdata") # load Rdata to global env

lng = c(-116.3529,-116.3726,-115.8352,-115.6001)
lat = c(33.77381,33.65167,33.48859,33.28386)
names = c("PD", "BDC", "DP","WI")

colorset6
mapviewOptions(fgb = FALSE)

# Create a data frame with the point data
points_df = data.frame(lng, lat, names)

# Convert the data frame to a spatial points data frame
points_sdf = st_as_sf(points_df,
                      coords = c("lng", "lat"), crs = 4326)

# Create color palette
pal<-c("#eb5e28","#390099","#ffbd00","#008000")

# Plot the points on a map
mapview(points_sdf, label = points_sdf$names,col.regions=pal,maxpixels = mapviewGetOption("mapview.maxpixels"),
        col=list("#eb5e28","#390099","#ffbd00","#008000"),legend=FALSE)

# save to object for file saving
ss.sites<-mapview(points_sdf, label = points_sdf$names,col.regions=pal,
                  col=list("#eb5e28","#390099","#ffbd00","#008000"),cex=8,legend=FALSE)


# NOTE: had to rearrange colors in legend in Photoshop, ran out of time to figure out in R
# html_fl = tempfile(fileext = ".html")
# png_fl = tempfile(fileext = ".png")

# save object to png
mapshot2(ss.sites, file = "figures/SaltonSeaStudySites3.png",
         vwidth=700,vheight=600,zoom=4,remove_controls=c("zoomControl", "layersControl", "homeButton"))
