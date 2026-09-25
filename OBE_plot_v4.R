library(viridis)
POsetDR.plot<-function(result, #output of DimensionalityReduction
                       N, #number of statistical units
                       k,#number of binary variables
                       gV_x=ceiling(k/2),#number of binary variables considered in building the grid on the x-axis
                       gV_y=k-ceiling(k/2),#number of binary variables considered in building the grid on the y-axis
                       Grid = TRUE, #if TRUE a grid based on gV variables is added on both axes
                       Grid_lty = c(1,5,2,4,3)[1:gV_x], #Line type of the grid segments associated to each variable in the grid
                       Grid_lwd = c(3,2.5,2,1.5,1)[1:gV_x], #Line width of the grid segments associated to each variable in the grid
                       Grid_col = "blue", #color of the grid
                       Profile = TRUE, #if TRUE a legend is added on both axes in order to recognize the profiles associated to each square in the grid
                       Profiles_col = "black", #color of the profile legend
                       varNames = NULL, #string /expression vector with the variables' names
                       Profiles_expand=c(1,1), # two dimension vector. The firts element regards x axis, the second y axis.
                       #if greater (lower) than one, the space between the lines in the profile legend is expanded (reduced).
                       title = "Optimal profile map", #main title of the plot
                       type=c("point", "jitter_square", "jitter_circle", "bubbles", "bubbles_and_approximation"), #type of plot
                       expandX = c(0,0), #add extra space to the default x-axis: c[1]=space added to the left, c[2]=space added to the right
                       expandY = c(0,0), #add extra space to the default y-axis: c[1]=space added to the left, c[2]=space added to the right
                       expand_radius=1, #add enalrge the radius of bubbles
                       points_col="black", #when type = "point", "jitter_square", "jitter_circle", it indicates the color of the points drawn, when
                       # type= "bubbles", "bubbles_and_approximation" it indicates the color of the border of the bubbles
                       center_col= "red", #when type = "jitter_square", "jitter_circle","bubbles", "bubbles_and_approximation"
                       #it indicates the color of the points drawn in the center of the clouds of points/bubbles
                       optimize_radius=FALSE,
                       bubbles_col = "gray", #when type ="bubbles" it indicates the bubbles' color
                       points_pch = 16, #when type = "point", "jitter_square", "jitter_circle", it indicates the type of the points drawn
                       center_pch = 16, #when type = "jitter_square", "jitter_circle", "bubbles", "bubbles_and_approximation",
                       #it indicates the type of the points drawn at the center of the clouds of points/bubbles
                       points_cex = 0.5, #when type = "point", "jitter_square", "jitter_circle", it indicates the dimension of the points drawn
                       center_cex = 0.5, #when type = "jitter_square", "jitter_circle", "bubbles", "bubbles_and_approximation",
                       #it indicates the dimension of the points drawn at the center of the clouds of points/bubbles
                       jitter_amount = 1, #when type = "jitter_square", "jitter_circle" it is proportional to the dispersion of the clouds of points
                       #around the centers
                       rampPalette = "viridis",
                       bubble_sizes = NULL,
                       cex.profile=1){
  #default variables' names: V[i] is the variable on the i-th column on the original dataset
  if(length(result)==3){
    names(result)[c(2,3)]<-c("bestVariablesPriority","bestRepresentation")}
  if(is.null(varNames)){
    varNames<-vector()
    for(i in 1:k){
      varNames<-append(varNames, parse(text=(paste0("V[",i,"]"))))
    }
  }
  #permuting variables names according to variables priority
  varNames<-varNames[result$bestVariablesPriority]
  
  
  #DIMENSIONS OF PLOT AREA
  #x and y lower bound: it depends on the presence/non presence of the profiles' legend in the plot
  ylow<-(-2^(k-gV_y)*gV_y)*(Profile == TRUE)+(Profile != TRUE)-expandY[1]
  xlow<-(-2^(k-gV_x)*gV_x)*(Profile == TRUE)+(Profile != TRUE)-expandX[1]
  #x and y upper bound: it depends on the presence/non presence of the profiles' legend in the plot
  #NB: type "bubbles and approximation" is a special case for xup
  yup<- (2^k+2^(k-gV_y)*3/2)*(Profile == TRUE)+2^k*(Profile != TRUE)+expandY[2]
  xup<- 2^k*(type != "bubbles_and_approximation")+(2^k+2^(k-gV_x)*3)*(type == "bubbles_and_approximation")+expandX[2]
  
  #DRAW THE PLOT AREA
  #x11()
  plot(x=NULL,xlim=c(xlow,xup), ylim=c(ylow,yup),xlab=" ", ylab=" ", main=title, yaxt="n",xaxt="n",bty="n")
  
  #ADDING THE POINTS (with optional jittering or bubbles) TO THE PLOT AREA
  xcoord<-result$bestRepresentation[,2]+0.5*(result$bestRepresentation[,2]==1)-0.5*(result$bestRepresentation[,2]==2^k)
  ycoord<-result$bestRepresentation[,3]+0.5*(result$bestRepresentation[,3]==1)-0.5*(result$bestRepresentation[,3]==2^k)
  if(type=="point"){
    points(xcoord, ycoord, pch=points_pch, cex=points_cex, col=points_col)
  }
  if(type=="jitter_square"){
    npoints<-round(N*result$bestRepresentation[,4]/sum(result$bestRepresentation[,4]))
    #npoints<-npoints+(npoints==0)
    x_unit<-rep(xcoord,npoints)
    y_unit<-rep(ycoord,npoints)
    Coordin1<-jitter(x_unit,amount=jitter_amount)*(rep(npoints,npoints)>1)+x_unit*(rep(npoints,npoints)==1)
    Coordin2<-jitter(y_unit,amount=jitter_amount)*(rep(npoints,npoints)>1)+y_unit*(rep(npoints,npoints)==1)
    points(Coordin1, Coordin2, pch=points_pch,  cex=points_cex, col=points_col)
    points(xcoord, ycoord, pch=center_pch,  cex=center_cex, col=center_col)
  }
  if(type=="jitter_circle"){
    npoints<-round(N*result$bestRepresentation[,4]/sum(result$bestRepresentation[,4]))
    #npoints<-npoints+(npoints==0)
    x_unit<-rep(xcoord,npoints)
    y_unit<-rep(ycoord,npoints)
    angolo<-runif(sum(npoints),max=2*pi)
    radius<-runif(sum(npoints), max=jitter_amount)
    Coordin1<-(x_unit+cos(angolo)*radius)*(rep(npoints,npoints)>1)+x_unit*(rep(npoints,npoints)==1)
    Coordin2<-(y_unit+sin(angolo)*radius)*(rep(npoints,npoints)>1)+y_unit*(rep(npoints,npoints)==1)
    points(Coordin1, Coordin2, pch=points_pch,  cex=points_cex, col=points_col)
    points(xcoord, ycoord, pch=center_pch,  cex=center_cex, col=center_col)
  }
  if(type=="bubbles"){
    if(optimize_radius==TRUE){
    DIStz<-(outer(xcoord,xcoord, function(x,y) (x-y)^2)+outer(ycoord,ycoord, function(x,y) (x-y)^2))^0.5
    diag(DIStz)<-max(DIStz)
    rad<-min(DIStz)}
    else{
      rad<-1
    }
    raggi<-(result$bestRepresentation[,4]/pi)^0.5
    raggi<-raggi*rad*expand_radius#raggi/max(raggi)*rad*expand_radius
    symbols(xcoord, ycoord, circles = raggi, inches = FALSE, add = TRUE, fg = points_col, bg = bubbles_col)
    points(xcoord, ycoord, pch=center_pch, cex=center_cex,  col=center_col)
  }
  if(type=="bubbles_and_approximation"){
    error<-result$bestRepresentation[,5]
    if(optimize_radius==TRUE){
      DIStz<-(outer(xcoord,xcoord, function(x,y) (x-y)^2)+outer(ycoord,ycoord, function(x,y) (x-y)^2))^0.5
      diag(DIStz)<-max(DIStz)
      rad<-min(DIStz)}
    else{
      rad<-1
    }
    raggi<-(result$bestRepresentation[,4]/pi)^0.5
    raggi<-raggi*rad*expand_radius#raggi/max(raggi)*rad*expand_radius
    
    n_col<-200
    colfunc <- viridis(n_col, option=rampPalette)
    scala_colori <- viridis(n_col, option=rampPalette)
    colori_ordine <- scala_colori[ceiling(error / max(error) * n_col) + (error == 0)]
    symbols(xcoord, ycoord,circles = raggi,inches = FALSE,add = TRUE,fg = points_col,bg = colori_ordine)
    points( xcoord, ycoord,  pch = center_pch,  cex = center_cex,  col = center_col)
    
    legend_image <- as.raster(matrix(viridis(n_col,option=rampPalette), ncol = 1))[n_col:1]
    text(x=2^k+2^(k-gV_x)*2.5, y = seq(0,2^k,l=5), labels = round(seq(0,max(error),l=5),2))
    rasterImage(legend_image, 2^k+2^(k-gV_x)*0.75+expandX[2], 1, 2^k+2^(k-gV_x)*1.75+expandX[2],2^k)
    }
  #ADDING THE (optional) GRID AND (optional) PROFILES' LEGEND
  #external borders of the grid
  segments(c(0,0,2^k,2^k), c(0,0,0,2^k),c(0,2^k,2^k,0), c(2^k,0,2^k,2^k),lwd=Grid_lwd[1], lty=Grid_lty[1], col= Grid_col)
  #grid on the x-axis
  for(i in 1:gV_x){
    v<-2^(k-i)*seq(1,2^i,2)+0.5
    l<-rep(0,length(v))
    L<-rep(2^k,length(v))
    if(Profile==TRUE){
      l2<-rep(-2^(k-gV_x)*Profiles_expand[1],length(v))
      L2<-rep(-2^(k-gV_x)*(gV_x+1)*Profiles_expand[1],length(v))
      segments(v,l,v,L2,lwd=1, lty=3, col=Profiles_col)
      axis(1, pos=-2^(k-gV_x)*(i-1)*Profiles_expand[1], lwd=0 ,at=c(1-2^(k-gV_x)/k,seq(2^(k-gV_x),2^k-1,2^(k-gV_x))+0.5,2^k)-2^(k-gV_x)/2,
           labels=c(varNames[gV_x-i+1],rep(c(0,1), times=2^(gV_x-i), each=2^(i-1))), col.axis=Profiles_col, cex.axis=cex.profile, gap.axis=0.01)
    }
    if(Grid==TRUE){
      segments(v,l,v,L,lwd=Grid_lwd[i], lty=Grid_lty[i], col=Grid_col)
    }
  }
  #grid on the y-axis
  for(i in 1:gV_y){
    h<-2^(k-i)*seq(1,2^i,2)+0.5
    l<-rep(0,length(h))
    L<-rep(2^k,length(h))
    if(Profile==TRUE){
      l2<-rep(-2^(k-gV_y)*Profiles_expand[2],length(h))
      L2<-rep(-2^(k-gV_y)*(gV_y+1)*Profiles_expand[2],length(h))
      segments(l,h,L2,h,lwd=1, lty=3, col=Profiles_col)
      axis(2, pos=-2^(k-gV_y)*(i-1)*Profiles_expand[2], lwd=0 ,at=c(seq(2^(k-gV_y),2^k-1,2^(k-gV_y))+0.5,2^k,2^k+2^(k-gV_y)*(1+1/k))-2^(k-gV_y)/2,
           labels=c(rep(c(0,1), times=2^(gV_y-i), each=2^(i-1)),varNames[gV_x+i]),col.axis=Profiles_col, cex.axis=cex.profile, gap.axis=0.01)
    }
    if(Grid==TRUE){
      segments(l,h,L,h,lwd=Grid_lwd[i], lty=Grid_lty[i], col=Grid_col)
    }
  }
  points(xcoord, ycoord, pch=points_pch, cex=points_cex, col=points_col)
}

profiles_freq_dist_weighted<-function(matrix, weights){
  #matrix è il dataframe che contiene la distribuzione di unità dei profili
  #osservati
  #weights è il vettore dei pesi di campionamento
  k<-ncol(matrix)
  #vettore dei profili decimalizzati
  dec_profiles<-as.matrix(matrix)%*%(2^((k:1)-1))
  observed_profiles_list<-sort(unique(dec_profiles))
  freq<-freq_W<-vector(length=length(observed_profiles_list))
  for(i in 1:length(observed_profiles_list)){
    freq_W[i]<-sum(weights*(dec_profiles==observed_profiles_list[i]))
    freq[i]<-sum(dec_profiles==observed_profiles_list[i])
  }
  obs_Profiles<-t(sapply(observed_profiles_list ,function(x){ as.integer(intToBits(x))})[k:1,])
  return(list(obs_Profiles,freq,freq_W))
}
