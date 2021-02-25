setwd('/Users/Avril/Desktop/all_ch_3_dirs/mapping/')

library(rgdal)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(rgeos)
library(ggplot2)
library(dplyr)
library(ggsn)
library(raster)
library(smoothr)

## create US-Canada border
can <- ne_countries(country='Canada', returnclass = 'sf', scale=10)
us <- ne_countries(country='United States of America', returnclass = 'sf', scale=10)
us.can <- st_intersection(us, can)
smooth.us.can <- smoothr::smooth(us.can, method='densify')
border = st_cast(smooth.us.can, to='POLYGON')

## get list of states and provinces to outline
states <- ne_states(country=c('United States of America', 'Canada'), returnclass = 'sf')
state.list <- c('Maine','Massachusetts','Vermont','New York','New Hampshire','Connecticut','Rhode Island',
                'Prince Edward Island','Ontario','Nova Scotia','Québec','New Brunswick','Pennsylvania',
                'New Jersey')
states <- states[which(states$name %in% state.list),]

## get river and lake data
rivs <- readOGR('/Users/Avril/Desktop/all_ch_3_dirs/mapping/Lakes_and_Rivers_Shapefile/NA_Lakes_and_Rivers/data/hydrography_l_rivers_v2.shp', stringsAsFactors=FALSE)
rivs <- st_as_sf(rivs)
rivs <- rivs[which(rivs$TYPE==17),]
LH <- rivs[which(rivs$NAMEEN=='LaHave River'),]

lakes <- readOGR('/Users/Avril/Desktop/all_ch_3_dirs/mapping/Lakes_and_Rivers_Shapefile/NA_Lakes_and_Rivers/data/hydrography_p_lakes_v2.shp', stringsAsFactors=FALSE)
lakes <- st_as_sf(lakes)
lakes <- lakes[which(lakes$TYPE==16),]
SE <- lakes[which(lakes$NAMEEN=="Sebago Lake"),]
LC <- lakes[which(lakes$NAMEEN=='Lake Champlain'),]

## plot map
al <- 1
line.w <- 0.5
green <- "gray95" 
blue <- '#DFF7FF' 
pdf('/Users/Avril/Desktop/map.pdf', width=6, height=4)
ggplot(data = states) +
  theme_classic() +
  theme(panel.background = element_rect(fill=alpha(blue, al)),
        axis.text=element_text(size=10, colour='black'),
        axis.title=element_blank(),
        panel.border=element_rect(fill=NA)) +
  geom_sf(fill=green, col='black', lwd=0.25) + 
  geom_sf(data = border, col='black', lwd=0.75) +
  ## add study waterbodies
  geom_sf(data = SE, fill=alpha(blue, al), col='#992913', lwd=line.w) + 
  geom_sf(data = LC, fill=alpha(blue, al), col='#1571a9', lwd=line.w) +
  geom_sf(data = LH, fill=NA, col='#D6A929', lwd=line.w) +
  ## add scalebar
  ggsn::scalebar(x.min=-74.5, x.max=-63.5, y.min=42.3, y.max=46.5,
           dist=50, dist_unit='km', st.bottom=TRUE, st.color='black',
           transform=TRUE, location='bottomright', st.size=2.5, st.dist=.03) +
  coord_sf(xlim = c(-74.5, -63), ylim = c(42, 46.5), expand = FALSE)
dev.off()

##### Lake Champlain map with greater extent for talks
pdf('/Users/Avril/Desktop/talk_map.pdf', width=5, height=5)
ggplot(data = states) +
  theme_classic() +
  theme(panel.background = element_rect(fill=alpha(blue, al)),
        axis.text=element_text(size=10, colour='black'),
        axis.title=element_blank(),
        panel.border=element_rect(fill=NA)) +
  geom_sf(fill=green, col='black', lwd=0.25) + ## other greens: springgreen4  #B5D383
  geom_sf(data = border, col='black', lwd=0.5) +
  ## add study waterbodies
  geom_sf(data = LC, fill=alpha(blue, al), col='#1571a9', lwd=line.w) +
  ## add all rivers
  # geom_sf(data = rivs, fill=NA, col='cadetblue1') +
  ## add scalebar
  # ggsn::scalebar(x.min=-80, x.max=-66.5, y.min=40, y.max=48,
  #                dist=50, dist_unit='km', st.bottom=TRUE, st.color='black',
  #                transform=TRUE, location='bottomright', st.size=2.5, st.dist=.03) +
  coord_sf(xlim = c(-80, -66.5), ylim = c(40, 48), expand = FALSE)
dev.off()
