#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Modeling Pb in soil surrounding Notre-Dame areas 
# copyright Yuling Yao and  Lex van Geen
# May 2020
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setwd("./ND") # set your dir and put the stan file that dir
soil.data=read.csv("soil2.csv")
x1=soil.data$Longitude
x2=soil.data$Latitude
y=soil.data$Soil.Pb..mg.kg.
plume=as.integer( soil.data$Relative.to.plume)==2
center=c(2.3496642, 48.8530245)# (Longitude , Latitude) of ND
type=as.integer( soil.data$Type)-1
type_word=levels(soil.data$Type)[-1]
library(RColorBrewer)
library(grDevices)
library(rstan)
options(mc.cores=8)
map2color<-function(x,pal,limits=NULL, breaks=NULL  ){
	if(is.null(limits)) limits=range(x)
	if(is.null(breaks ))  breaks=seq(limits[1],limits[2],length.out=length(pal)+1)
	pal[findInterval(x,  breaks , all.inside=TRUE)]
}
my.palette <- colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(20)  # generate color shades
col_brew=map2color( log10(y),my.palette, breaks = seq(1,4,length.out = 20))
library("geosphere")
n=length(y)
d=a=a_radian2=c()   ##  distance and angle of all collected samples
for( i in 1:n){
	d[i]=	distGeo(c(center[1],  center[2]), c(x1[i], x2[i]))
	#a[i]=angleFromCoordinate(center[2], center[1], x2[i], x1[i] )
	a[i]= (bearing(c(center[1],  center[2]), c(x1[i], x2[i]))+ 360) %%360
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Fit a Gaussian process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
d_unit=d/1000   ## convert m to KM
y_log=log10(y) ## use log (base 10) in the actual model 
a_radian=a/180*pi ## convert degree from (0,360) to (0, 2pi)
circle= (as.integer( soil.data$Original.circles)==3)



n_0= 29     ## number of grids in each dimension.
grid_d=seq(0.1, 1.5,length.out = n_0)
grid_theta=seq(0, 2*pi,length.out = n_0+1)[-c(n_0+1)]
delta_d=grid_d[2]-grid_d[1]  ## spacing
delta_theta=grid_theta[2]-grid_theta[1]
temp=grid_1Dto2D=matrix(NA,  2, n_0^2)
for (i in c(1:n_0))
	for (j in c(1:n_0)){
		temp[,(i-1)*n_0+j]= c(grid_d[i], grid_theta[j])
		grid_1Dto2D[,(i-1)*n_0+j]= c(i,j)
	}
d_test=temp[1,]
theta_test=temp[2,]
m2= stan_model("type.stan")
y_forth=y^(1/4)
stan_fit=sampling (m2, data=list(N=length(y_forth), y=y_forth, 
																 d=d_unit, theta=a_radian, 
																 N_test=length(d_test), 
																 d_test=d_test, theta_test=theta_test,
																 type=type),
									 iter=3000,chains=4)
print(stan_fit)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Posterior inference     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


fit_draw=rstan::extract(stan_fit)
f_predict_mean=colMeans(fit_draw$f_predict)
f_predict_sd=apply(fit_draw$f_predict, 2, sd)
range(f_predict_sd)
range(f_predict_mean)^4
range(log10(f_predict_mean^4))
col_brew_test=map2color(f_predict_mean^4 ,my.palette, breaks = seq(30, 650 ,length.out = 20))
col_brew_test_sd=map2color(f_predict_sd ,my.palette, breaks = seq(0.2, 0.8,length.out = 20))
col_brew_test_log=map2color( log10(f_predict_mean^4), my.palette, breaks = seq(log10(30), log10(650), length.out = 20))
f_draw=(fit_draw$f_predict)^4
y_draw=(fit_draw$y_predict)^4


range(a[plume])
  ## use the data label to determine the boundary of the plume : 260 < angle < 310
id_in_test= which(theta_test*180/pi >259.5 & theta_test*180/pi <310.5   ) 
id_in_test_grid= which(grid_theta*180/pi >259.5 & grid_theta*180/pi <310.5   ) 
which_id_is_grid_d=function(i){
	which(grid_1Dto2D[1,]== i)
}
which_id_smaller_grid_d=function(i){
	which(grid_1Dto2D[1,]<= i)
}

which_id_is_grid_theta=function(i){
	which(grid_1Dto2D[2,]== i)
}

S= dim(f_draw)[1]
average_across_distance=average_across_theta=average_across_theta_d_500=average_across_theta_d_1000=matrix(NA, S, n_0)
for (i in 1:n_0)
	average_across_distance[,i]=  rowMeans(f_draw [,which_id_is_grid_d(i)])

for (i in 1:n_0){
	average_across_theta[,i]=  rowMeans(f_draw [,which_id_is_grid_theta(i)])
	average_across_theta_d_500[,i]=  rowMeans(f_draw [,	which(grid_1Dto2D[2,]== i ) ][,5:9] )  
	# print( grid_d[5:9] )
	average_across_theta_d_1000[,i]=  rowMeans(f_draw [,	which(grid_1Dto2D[2,]== i)  ][,17:21]  ) 
	# print( grid_d[17:21] )																												
}
colMeans(average_across_theta_d_500[,id_in_test_grid])
mean(  rowMeans(average_across_theta_d_500[,-id_in_test_grid]))
quantile(  rowMeans(average_across_theta_d_500[,id_in_test_grid]), c(0.025,0.975))
apply(average_across_theta_d_500[,id_in_test_grid], 2, mean )
apply(average_across_theta_d_500[,id_in_test_grid], 2, quantile,  c(0.025,0.975))
average_across_plume=average_across_nplume=average_across_plume_diff=y_average_across_plume_diff=y_average_across_plume_diff_exp=average_across_plume_diff_exp=average_across_plume_exp=average_across_nplume_exp=y_average_across_plume_diff_anti_exp=average_across_plume_diff_anti_exp=average_across_plume_diff_exp_inside=y_average_across_plume_diff_exp_inside=average_across_plume_exp_inside=average_across_nplume_exp_inside=average_across_plume_y=average_across_nplume_y=matrix(NA, S, n_0)
for (i in 1:n_0){
	average_across_plume[,i]=  rowMeans((f_draw [,which_id_is_grid_d (i)[id_in_test_grid ]  ])) 
	average_across_nplume[,i]=  rowMeans((f_draw [,which_id_is_grid_d (i)[-id_in_test_grid ]  ])) 
	average_across_plume_y[,i]=  rowMeans((y_draw [,which_id_is_grid_d (i)[id_in_test_grid ]  ])) 
	average_across_nplume_y[,i]=  rowMeans((y_draw [,which_id_is_grid_d (i)[-id_in_test_grid ]  ])) 
	average_across_plume_diff[,i]=  rowMeans(f_draw [,which_id_is_grid_d (i)[id_in_test_grid ]  ]) -rowMeans(f_draw [,which_id_is_grid_d (i)[-id_in_test_grid]  ]) 
	average_across_plume_diff_exp_inside[,i]=  rowMeans(( f_draw [,intersect(id_in_test, which_id_smaller_grid_d (i))])  )   -rowMeans(( f_draw [,intersect( c(1:(n_0^2)) [-id_in_test], which_id_smaller_grid_d (i))])  ) 
	y_average_across_plume_diff[,i]=  rowMeans(y_draw [,which_id_is_grid_d (i)[id_in_test_grid ]  ]) -rowMeans(y_draw [,which_id_is_grid_d (i)[-id_in_test_grid]  ]) 
}
 excess_area=excess_area_y=excess_in_circle= excess_in_circle_y=weighted_plume=weighted_nplume=weighted_plume_y=weighted_nplume_y=matrix(NA, dim (average_across_plume_diff)[1],  n_0)
for( j in 1:n_0)
	for( s in 1:dim (average_across_plume_diff)[1] ){
		grid_d_s=grid_d[1:j]
		excess_area[s, j]= sum( average_across_plume_diff[s, 1:j] * grid_d_s   /sum(grid_d_s))
		excess_area_y[s, j]= sum( y_average_across_plume_diff[s,1:j] * grid_d_s /sum(grid_d_s))
		excess_in_circle[s,j]  =pi  *  (grid_d[j] * 1000)^2 *  6/60  * 0.01   *2* 1000 * excess_area[s, j]  /1e6
		excess_in_circle_y[s,j]  =pi  *  (grid_d[j] * 1000)^2 *  6/60  * 0.01   *2* 1000 * excess_area_y [s, j] /1e6
		weighted_plume[s,j]=sum( average_across_plume[s, 1:j] * grid_d_s   /sum(grid_d_s))
		weighted_nplume[s,j]=sum( average_across_nplume[s, 1:j] * grid_d_s   /sum(grid_d_s))
		weighted_plume_y[s,j]=sum( average_across_plume_y[s, 1:j] * grid_d_s   /sum(grid_d_s))
		weighted_nplume_y[s,j]=sum( average_across_nplume_y[s, 1:j] * grid_d_s   /sum(grid_d_s))
}
 
type_draw=extract(stan_fit, pars="type_eff")$type_eff



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Prediction in the ordinal coordinate %%%%%%%%%%%%%%%% 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

n_0=31
x1_grid=seq(-1.5, 1.5, length.out =n_0)
x2_grid=seq(-1.5, 1.5, length.out =n_0)
x1_data= d* sin(a/180*pi)
x2_data= d*cos(a/180*pi)

delta_x=x1_grid[2]-x1_grid[1]
temp=matrix(NA,  2, n_0^2)
for (i in c(1:n_0))
	for (j in c(1:n_0))
		temp[,(i-1)*n_0+j]= c(x1_grid[i], x2_grid[j])
x_grid=temp
x1_test=x_grid[1,]
x2_test=x_grid[2,]
d_test_new=sqrt(x1_test^2 +x2_test^2 )
theta_test_new=pi/2- atan2( x2_test, x1_test)
 

stan_fit=sampling (m2, data=list(N=length(y_forth), y=y_forth, 
																 d=d_unit, theta=a_radian, 
																 N_test=length(d_test_new), 
																 d_test=d_test_new, theta_test=theta_test_new,
																 type=type),
									 iter=3000,chains=4)
print(stan_fit)
 
fit_draw=extract(stan_fit)
f_predict_mean=apply(  (fit_draw$f_predict)^4 , 2, mean)
f_predict_mean_matrix=matrix(NA, n_0, n_0)
for (i in c(1:n_0))
	for (j in c(1:n_0))
		f_predict_mean_matrix[i,j]= f_predict_mean[(i-1)*n_0+j]
