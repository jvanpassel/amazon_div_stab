
## Load alpha stability: EVI TAC trend
evi<-brick('evi_01_19_masked.tif')
evi_tac_60m_slope<-raster('evi_01_19_lag1_60m.tif')
evi_tac_60m_slope_df<-read.csv('evi_01_19_lag1_60m_nonaslope.csv')[,-1]

## Add environmental variables ####
dem_amazon<-raster('SRTM_DEM_amazon.tif')
interann_var<-raster('betweenyearvariation_80_19.tif')
seasonality<-raster('seasonality_index_80_19.tif')
soil_s_0_100<-raster('soil_sand_0_100_amazon.tif')
soil_c_0_100<-raster('soil_clay_0_100_amazon.tif')
drought<-read.csv('droughtlegacy_01_19_nonaslope.csv')[,-1]
drought[is.na(drought)]<-0
evi_tac_60m_slope_df$elevation<-raster::extract(dem_amazon,evi_tac_60m_slope_df[,1:2])
evi_tac_60m_slope_df$interann_var<-raster::extract(interann_var,evi_tac_60m_slope_df[,1:2])
evi_tac_60m_slope_df$seasonality<-raster::extract(seasonality,evi_tac_60m_slope_df[,1:2])
evi_tac_60m_slope_df$s_0_100<-raster::extract(soil_s_0_100,evi_tac_60m_slope_df[,1:2])
evi_tac_60m_slope_df$c_0_100<-raster::extract(soil_c_0_100,evi_tac_60m_slope_df[,1:2])
evi_tac_60m_slope_df<-cbind(evi_tac_60m_slope_df,drought[,c(3,5,6,11)])

## Load Keil and Chase diversity data, cropped to Amazon ####
div_kc_1ha_amazon<-shapefile('pred_1haplots_amazon.shp')
div_kc_grids_amazon<-shapefile('pred_grids_amazon.shp')

## Add alpha, beta, gamma diversity per pixel ####
#change crs of coordinates and points to UTM (meters)
stab_df<-evi_tac_60m_slope_df
coordinates(stab_df)<-~x+y
crs(stab_df)<-crs(evi_tac_60m_slope)
stab_df_utm<-spTransform(stab_df,'+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs ')
div_kc_1ha_amazon_utm<-spTransform(div_kc_1ha_amazon,'+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs ')
#convert to sf 
stab_df_utm_sf<-st_as_sf(stab_df_utm)
div_kc_1ha_amazon_utm_sf<-st_as_sf(div_kc_1ha_amazon_utm)

#fill in mean value per buffer surrounding diversity point
df <- list()
for(i in 1:nrow(div_kc_1ha_amazon_utm_sf)){
  pointi<-div_kc_1ha_amazon_utm_sf[i,]
  #create buffer of 8500/14100/19700 km to get area of 225/625/1225 km²
  pointi_buffer<-st_buffer(pointi$geometry,8500)
  #see which points intersect with the buffer zone and extract coordinates
  int_buffer<-st_intersects(stab_df_utm_sf,pointi_buffer,sparse=F)
  tac_in_buffer<-stab_df_utm_sf[which(int_buffer),]
  #create dataframe with all env and drought var + diversity per drought
  if(!all(is.na(tac_in_buffer))){
    tac_in_buffer<-cbind(tac_in_buffer,as.data.frame(st_coordinates(tac_in_buffer)))
    tac_in_buffer<-as.data.frame(tac_in_buffer)[,-c(16)]
    coordinates(tac_in_buffer)<-~X+Y
    crs(tac_in_buffer)<-'+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs '
    tac_in_buffer_wgs<-spTransform(tac_in_buffer,crs(evi_tac_60m_slope))
    tac_in_buffer_wgs<-as.data.frame(tac_in_buffer_wgs)
    #create matrix with equal amount of rows with diversity data
    pointi_div<-as.numeric(pointi[,c(15,20,22,17)])[-5]
    pointi_div_matrix<-rbind(pointi_div)[rep(1,nrow(tac_in_buffer_wgs)),]
    if(nrow(tac_in_buffer_wgs)==1){
      tac_in_buffer_wgs<-cbind(tac_in_buffer_wgs,as.data.frame(matrix(pointi_div_matrix,nrow=1)))
    } else {
      tac_in_buffer_wgs<-cbind(tac_in_buffer_wgs,pointi_div_matrix)
    }
    tac_in_buffer_wgs$diversity_id<-i
    colnames(tac_in_buffer_wgs)[14:17,19:20]<-c("x","y","alphadiv","betadiv_loc",
                                          "gammadiv","diversity_id")
    df[[i]] <- tac_in_buffer_wgs
  }
  print(i)
}
df_list<-df
df_list<-do.call('rbind',df_list)
df_list<-as.data.frame(df_list)
colnames(df_list)[1]<-'slope'
write.csv(df_list,'evi_01_19_60m_tac_slope_withenv_buffer_indpixels_225km2_kc_div_wgs.csv')[,-1]

## Calc gamma stability from mean EVI time series per grid and then TAC slope ####
dates<-seq(as.Date('2001-01-01'),as.Date('2019-12-01'),by='month')
mo <- as.numeric(format(dates[1], "%m"))
yr <- as.numeric(format(dates[1], "%Y"))
t<-c(1:169)
grid_evi_mean<-as.data.frame(matrix(ncol=228,nrow=length(div_kc_grids_amazon)))
grid_tac_mean<-as.data.frame(matrix(ncol=169,nrow=length(div_kc_grids_amazon)))
grid_tac_slope<-as.data.frame(matrix(ncol=1,nrow=length(div_kc_grids_amazon)))

for(i in 1:length(div_kc_grids_amazon)){
  grid_i<-div_kc_grids_amazon[i,]
  grid_i_crs<-spTransform(grid_i,crs(evi))
  #crop EVI to grid cell
  evi_i<-crop(evi,grid_i_crs)
  evi_i<-mask(evi_i,grid_i_crs)
  #calc mean EVI per grid cell
  evi_i_df<-as.data.frame(evi_i)
  evi_i_mean<-colMeans(evi_i_df,na.rm=T)
  evi_i_mean[evi_i_mean<0]<-NA
  grid_evi_mean[i,]<-evi_i_mean
  #STL decomposition
  if(!all(is.na(evi_i_mean))){
    evi_i_mean<-na_interpolation(evi_i_mean,option='linear',maxgap=1)
    evi_i_ts<-ts(evi_i_mean, start = c(yr, mo), freq = 12)
    evi_i_stl<-stl(evi_i_ts,s.window = 13,t.window = 19,l.window = 13)
    #calc TAC from remainder
    evi_i_remain<-evi_i_stl$time.series[,3]
    evi_i_tac <- rollapply(evi_i_remain, width = 60,
                           FUN = function(z) acf(z,na.action=na.pass,lag.max= 1,plot=FALSE)$acf[2],
                           by.column = FALSE, align = "right")
    grid_tac_mean[i,]<-evi_i_tac
    #calc TAC slope
    evi_i_lm_slope<-lm(evi_i_tac~t)
    tac_slope_i<-evi_i_lm_slope$coefficients[[2]]
    grid_tac_slope[i,]<-tac_slope_i
  }
  print(i)
}
write.csv(grid_evi_mean,'evi_01_19_pergrid_meanevi.csv')
write.csv(grid_tac_mean,'evi_01_19_pergrid_meantac.csv')
write.csv(grid_tac_slope,'evi_01_19_pergrid_meantacslope.csv')

#link mean TAC slope values to grids
grid_tac_slope<-read.csv('evi_01_19_pergrid_meantacslope.csv')[,-1]
div_kc_grids_amazon$meantacslope<-grid_tac_slope
writeOGR(div_kc_grids_amazon,'pred_grids_amazon_meantacslope.shp',
         'pred_grids_amazon_meantacslope',driver='ESRI Shapefile',overwrite_layer = T)

## Calc beta stability on local scale ####
ind_tac<-read.csv('evi_01_19_60m_tac_slope_withenv_buffer_indpixels_225km2_kc_div_wgs.csv')[,-1]
coordinates(ind_tac)<-~x+y
crs(ind_tac)<-crs(evi_tac_60m_slope)
ind_tac<-st_as_sf(ind_tac)
div_kc_grids_amazon_gammastab<-shapefile('./keil and chase, 2019/pred_grids_amazon_meantacslope.shp')
div_kc_grids_amazon_gammastab<-spTransform(div_kc_grids_amazon_gammastab,crs(evi_tac_60m_slope))
div_kc_grids_amazon_gammastab_sf<-st_as_sf(div_kc_grids_amazon_gammastab)
ind_tac_inter<-st_intersection(ind_tac,div_kc_grids_amazon_gammastab_sf)
ind_tac_inter<-ind_tac_inter[,c(1:19,43)]
ind_tac_inter<-cbind(ind_tac_inter,as.data.frame(st_coordinates(ind_tac_inter)))
ind_tac_inter<-as.data.frame(ind_tac_inter)[,-c(23)]
colnames(ind_tac_inter)[19:22]<-c('region_id','gammastab','x','y')
ind_tac_inter$betastab_local<-(ind_tac_inter$gammastab+0.1)/(ind_tac_inter$slope+0.1)
write.csv(ind_tac_inter,'evi_01_19_60m_tac_slope_withenv_buffer_indpixels_225km2_kc_div_wgs_a_b_g.csv')

## Piecewise SEM ####
tac_ind_div<-read.csv('evi_01_19_60m_tac_slope_withenv_buffer_indpixels_225km2_kc_div_wgs_a_b_g.csv')
tac_ind_div<-tac_ind_div[complete.cases(tac_ind_div),]
tac_ind_div$x<-jitter(tac_ind_div$x)
tac_ind_div$y<-jitter(tac_ind_div$y)
tac_ind_div$diversity_id <- as.factor(tac_ind_div$diversity_id)
tac_ind_div$region_id <- as.factor(tac_ind_div$region_id)
tac_ind_div$slope <- (-tac_ind_div$slope)
tac_ind_div$gammastab <- (-tac_ind_div$gammastab)
tac_ind_div$betastab_local_div <- (tac_ind_div$gammastab + 0.1)/(tac_ind_div$slope + 0.1)
tac_ind_div$gammadiv_int<-as.integer(tac_ind_div$gammadiv)

modelList_all <- psem(
  lme(slope ~ alphadiv + elevation + s_0_100 + c_0_100 + d_amount_01_19 + 
        d_avdur_01_19 + d_avint_01_19, 
      random = ~1|diversity_id,correlation = corExp(form=~x+y),data = tac_ind_div),
  lme(betastab_local_div ~ betadiv_loc + seasonality + elevation + s_0_100 + c_0_100 + 
        d_amount_01_19 + d_avdur_01_19 + d_avint_01_19, 
      random = ~1|diversity_id/region_id,correlation = corExp(form=~x+y),data = tac_ind_div),
  lme(gammastab ~ slope + betastab_local_div + gammadiv_int, 
      random = ~1|diversity_id,correlation = corExp(form=~x+y),data = tac_ind_div),
  glmmPQL(gammadiv_int ~ alphadiv + betadiv_loc,family = poisson(link='log'),
          random = ~1|region_id,data = tac_ind_div),
  slope%~~%betastab_local_div,
  alphadiv%~~%betadiv_loc,
  alphadiv%~~%betastab_local_div,
  slope%~~%betadiv_loc,
  alphadiv%~~%seasonality,
  gammadiv_int%~~%elevation,
  gammadiv_int%~~%c_0_100,
  gammadiv_int%~~%seasonality,
  gammadiv_int%~~%d_amount_01_19,
  data = tac_ind_div
)
coefs(modelList_all)
summary(modelList_all) 
rsquared(modelList_all, method = NULL)
fisherC(modelList_all)
dSep(modelList_all)

