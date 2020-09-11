## weight TL rivers by HydroAtlas watershed area
## Sept 11, 2020
## C Kim

library('sf')
library('here')
library('tidyverse')

# read in sed plume Timor tif
# https://cfss.uchicago.edu/notes/simple-features/
pp <- here("Data/shape/NTimor_pour_points.shp") %>% 
   st_read()
   
pp
