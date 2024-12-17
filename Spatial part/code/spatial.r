


project_points = function(points, slope, intercept){
    apply(points, 1, function(a){
        x0 = a[[1]]
	y0 = a[[2]]
	x1 = (x0 + slope*y0 - slope*intercept)/(1 + slope**2)
	y1 = (slope*x0 + slope**2*y0 + intercept)/(1 + slope**2)
	d1 = sqrt((x1 - 0)**2 + (y1 - intercept)**2)
	return(c(x0=x0, y0=y0, x1=x1, y1=y1, d1=d1))
    })
}


fast_moran_i = function(x, w){
    # x = numeric vector
    # w = weights (can be sparse)
    x = x - mean(x)
    length(x)/sum(w)*sum(x * t(x * w))/sum(x**2)
}
