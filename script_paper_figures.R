library(raster)
library(sf)
library(scico)
library(khroma)
library(ggplot2)
library(zoo)
library(ggfortify)

evi<-brick('evi_01_19_masked.tif')
evi_tac_60m_slope<-raster('evi_01_19_lag1_60m.tif')
evi_tac_60m_slope_df<-read.csv('evi_01_19_lag1_60m_nonaslope.csv')[,-1]

tac_ind_div_225<-read.csv('evi_01_19_60m_tac_slope_withenv_buffer_indpixels_225km2_kc_div_wgs_a_b_g.csv')[,-1]
div_kc_grids_amazon_gammastab<-shapefile('./keil and chase, 2019/pred_grids_amazon_betadiv_reg_meantacslope.shp')

## Plot alpha-alpha and beta-beta relationships ####
tac_ind_div<-tac_ind_div_225
tac_ind_div<-tac_ind_div[complete.cases(tac_ind_div),]
tac_ind_div$x<-jitter(tac_ind_div$x)
tac_ind_div$y<-jitter(tac_ind_div$y)
tac_ind_div$diversity_id <- as.factor(tac_ind_div$diversity_id)
tac_ind_div$region_id <- as.factor(tac_ind_div$region_id)
tac_ind_div$slope<-(-tac_ind_div$slope)
tac_ind_div$gammastab<-(-tac_ind_div$gammastab)
tac_ind_div$betastab_local_div<-(tac_ind_div$gammastab+0.1)/(tac_ind_div$slope+0.1)
summary(lme(slope~alphadiv,data = tac_ind_div,method='ML',
            random = ~1|diversity_id, correlation = corExp(form=~x+y)))
summary(lme(betastab_local_div~betadiv_loc,data = tac_ind_div,method='ML',
            random = ~1|diversity_id/region_id, correlation = corExp(form=~x+y)))

vikO<-colour('vikO')
buda<-colour('buda')

#plot alpha-alpha
stderror_int<-2.431098e-04
val_int<-c(-2.681583e-04+stderror_int,-2.681583e-04-stderror_int)
stderror_slope<-1.968330e-06
val_slope<-c(6.080540e-06+stderror_slope,6.080540e-06-stderror_slope)

p <- ggplot(tac_ind_div, aes(alphadiv, slope)) + 
  geom_point(alpha=0.2) + 
  geom_abline(intercept=-2.681583e-04,slope=6.055510e-06,col='black') +
  ylim(c(-0.005,0.005))+
  theme_bw() +
  theme(text=element_text(size=17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab('Alpha diversity') +  ylab('Alpha stability')
p.x <- layer_scales(p)$x$get_limits()
p.y <- layer_scales(p)$y$get_limits()
df1 <- data.frame(x = rep(c(p.x[1] - (p.x[2] - p.x[1]),
                            p.x[2] + (p.x[2] - p.x[1])), each = 2),
                  intcpt = c(val_int, rev(val_int))
) %>%  mutate(y = intcpt + 6.055510e-06 * x)
df2 <- data.frame(x = rep(c(p.x[1] - (p.x[2] - p.x[1]),
                            p.x[2] + (p.x[2] - p.x[1])), each = 2),
                  slope = c(val_slope, rev(val_slope))
) %>%  mutate(y = -2.681583e-04 + slope * x)
p + 
  annotate(geom = "polygon",x = df1$x,y = df1$y,fill='lightgrey') +
  annotate(geom = "polygon",x = df2$x,y = df2$y,fill='lightgrey') +
  coord_cartesian(xlim = p.x, ylim = p.y)

ggplot(tac_ind_div, aes(alphadiv, slope)) + 
  annotate(geom = "polygon",x = df1$x,y = df1$y,fill='lightgrey') +
  annotate(geom = "polygon",x = df2$x,y = df2$y,fill='lightgrey') +
  coord_cartesian(xlim = p.x, ylim = p.y) +
  geom_point(alpha=0.2) + 
  geom_abline(intercept=-2.681583e-04,slope=6.055510e-06,col='black') +
  ylim(c(-0.005,0.005))+
  theme_bw() +
  theme(text=element_text(size=17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab('Alpha diversity') +  ylab('Alpha stability')

#combine plots
ggplot(tac_ind_div, aes(x = alphadiv, fill=..x..)) +
  geom_histogram(show.legend = F) +
  scale_fill_scico(palette = 'buda', direction = -1,
                   na.value = 'transparent') +
  theme_bw() + 
  theme(text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(tac_ind_div, aes(y = slope, fill=..y..)) +
  geom_histogram(show.legend = F) +
  scale_fill_scico(palette = 'vikO', direction = -1,
                   name='Alpha\nstability',
                   na.value = 'transparent',
                   limits=c(-0.005,0.005),midpoint=0) +
  theme_bw() + 
  theme(text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#plot beta-beta
stderror_int<-0.003215990 
val_int<-c(0.9976059+stderror_int,0.9976059-stderror_int)
stderror_slope<-0.000196859  
val_slope<-c(0.0008213+stderror_slope,0.0008213-stderror_slope)

p <- ggplot(tac_ind_div, aes(betadiv_loc, betastab_local_div)) + 
  geom_point(alpha=0.2) + 
  geom_abline(intercept=0.9976105,slope=0.0008200,col='black') +
  theme_bw() +
  theme(text=element_text(size=17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab('Beta diversity') +  ylab('Spatial asynchrony')
p.x <- layer_scales(p)$x$get_limits()
p.y <- layer_scales(p)$y$get_limits()
df1 <- data.frame(x = rep(c(p.x[1] - (p.x[2] - p.x[1]),
                            p.x[2] + (p.x[2] - p.x[1])), each = 2),intcpt = c(val_int, rev(val_int))
) %>%  mutate(y = intcpt + 0.0008200 * x)
df2 <- data.frame(x = rep(c(p.x[1] - (p.x[2] - p.x[1]),
                            p.x[2] + (p.x[2] - p.x[1])), each = 2),
                  slope = c(val_slope, rev(val_slope))
) %>%  mutate(y = 0.9976059 + slope * x)

p + 
  annotate(geom = "polygon",x = df1$x,y = df1$y,fill='lightgrey') +
  annotate(geom = "polygon",x = df2$x,y = df2$y,fill='lightgrey') +
  coord_cartesian(xlim = p.x, ylim = p.y)

#combine plots
ggplot(tac_ind_div, aes(betadiv_loc, betastab_local_div)) + 
  annotate(geom = "polygon",x = df1$x,y = df1$y,fill='lightgrey') +
  annotate(geom = "polygon",x = df2$x,y = df2$y,fill='lightgrey') +
  coord_cartesian(xlim = p.x, ylim = p.y) +
  geom_point(alpha=0.2) + 
  geom_abline(intercept=0.9976105,slope=0.0008200,col='black') +
  theme_bw() +
  theme(text=element_text(size=17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab('Beta diversity') +  ylab('Spatial asynchrony')

ggplot(tac_ind_div, aes(x = betadiv_loc, fill=..x..)) +
  geom_histogram(show.legend = F) +
  scale_fill_scico(palette = 'buda', direction = -1,
                   na.value = 'transparent') +
  theme_bw() + 
  theme(text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(tac_ind_div, aes(y = betastab_local_div, fill=..y..)) +
  geom_histogram(show.legend = F) +
  scale_fill_scico(palette = 'vikO', direction = -1,
                   name='Spatial',
                   na.value = 'transparent',
                   midpoint=1) +
  theme_bw() + 
  theme(text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Plot stability and diversity across spatial scales ####
tac_ind_div<-tac_ind_div_225
tac_ind_div<-tac_ind_div[complete.cases(tac_ind_div),]
tac_ind_div$slope<-(-tac_ind_div$slope)
tac_ind_div$gammastab<-(-tac_ind_div$gammastab)
tac_ind_div$betastab_local_div<-(tac_ind_div$gammastab + 0.1)/(tac_ind_div$slope + 0.1)

coordinates(tac_ind_div)<-~x+y
crs(tac_ind_div)<-crs(evi_tac_60m_slope)
tac_ind_div_sf<-st_as_sf(tac_ind_div)
div_kc_grids_amazon_sf<-st_as_sf(div_kc_grids_amazon_gammastab)

ggplot(tac_ind_div) +
  geom_histogram(aes(x=slope),color='black',fill='white',binwidth = 0.001)+
  theme_bw() +
  theme(text=element_text(size=30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab('Alpha stability') +
  ylab('Frequency')

vikO<-colour('vikO')
batlow<-colour('batlow')
buda<-colour('buda')

#alpha stab
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',
            fill='white',size=0.3) +
  geom_sf(data = div_kc_grids_amazon_sf, aes(geometry=geometry),color='black',linewidth=0.3) +
  geom_sf(data = tac_ind_div_sf, aes(color=slope),size=1,alpha=0.5) +
  scale_colour_scico(palette = 'vikO', direction = -1,
                     name='Alpha\nstability',
                     na.value = 'transparent',
                     limits=c(-0.005,0.005),midpoint=0) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name=expression(paste("Longitude [deg]")),
                     limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0))

#beta stab
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',
            fill='white',size=0.3) +
  geom_sf(data = div_kc_grids_amazon_sf, aes(geometry=geometry),color='black',linewidth=0.3) +
  geom_sf(data = tac_ind_div_sf, aes(color=betastab_local_div),size=1,alpha=0.5) +
  scale_colour_scico(palette = 'vikO', direction = -1,
                     name='Spatial',
                     na.value = 'transparent',
                     midpoint=1) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name=expression(paste("Longitude [deg]")),
                     limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0))

#gamma stab
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',
            fill='white',size=0.3) +
  geom_sf(data = div_kc_grids_amazon_sf, aes(geometry=geometry,fill=mentcslp)) +
  scale_fill_scico(palette = 'vikO', direction = -1,
                   name='Gamma\nstability',
                   na.value = 'transparent',
                   limits=c(-0.005,0.005),midpoint=0) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name=expression(paste("Longitude [deg]")),
                     limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0))

#alpha div
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',
            fill='white',size=0.3) +
  geom_sf(data = div_kc_grids_amazon_sf, aes(geometry=geometry),color='black',linewidth=0.3) +
  geom_sf(data = tac_ind_div_sf, aes(color=alphadiv),size=1) +
  scale_colour_scico(palette = 'buda', direction = -1,
                     name='Alpha',
                     na.value = 'transparent') +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name=expression(paste("Longitude [deg]")),
                     limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0))

#beta stab
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',
            fill='white',size=0.3) +
  geom_sf(data = div_kc_grids_amazon_sf, aes(geometry=geometry),color='black',linewidth=0.3) +
  geom_sf(data = tac_ind_div_sf, aes(color=betadiv_loc),size=1) +
  scale_colour_scico(palette = 'buda', direction = -1,
                     name='Spatial',
                     na.value = 'transparent') +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name=expression(paste("Longitude [deg]")),
                     limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0))

#gamma div
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',
            fill='white',size=0.3) +
  geom_sf(data = div_kc_grids_amazon_sf, aes(geometry=geometry,fill=S)) +
  scale_fill_scico(palette = 'buda', direction = -1,
                   name='Gamma',
                   na.value = 'transparent') +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name=expression(paste("Longitude [deg]")),
                     limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0))

