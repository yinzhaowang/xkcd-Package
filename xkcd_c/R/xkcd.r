
rxkcd <- function(n, sd=1){
  if(length(sd) == 1){
    sd = rep(sd, length(n))
  }
  if(length(sd) < length(n)){
    sd = append(sd, rep(1, length(length(n) - length(sd))))
  }
  
  out <- .Call("rxkcd", n = as.integer(n), sd = as.double(sd) , PACKAGE = "xkcd")
  return(out)
}

dxkcd <- function(x, sd = 1, log = FALSE, swap.end.points = FALSE){
  if(length(sd) == 1){
    sd = rep(sd, length(x))
  }
  if(length(sd) < length(x)){
    sd = append(sd, rep(1, length(length(x) - length(sd))))
  }
  
  
  out <- .Call("dxkcd", r_x = as.double(x), r_d = as.double(sd), r_log= as.integer(log), 
        r_swap_end_points= as.integer(swap.end.points), r_n = length(x), r_sd_n = length(sd) , PACKAGE = "xkcd")
  return(out)
}

pxkcd <- function(q, sd = 1, log.p = FALSE, swap.end.points = FALSE){
  if(length(sd) == 1){
    sd = rep(sd, length(q))
  }
  if(length(sd) < length(q)){
    sd = append(sd, rep(1, length(length(q) - length(sd))))
  }
  
  out <- .Call("pxkcd", r_x = as.double(q), r_d = as.double(sd), r_log= as.integer(log.p), r_swap_end_points= as.integer(swap.end.points)
      , r_n = length(q), r_sd_n = length(sd) , PACKAGE = "xkcd")
  return(out)
}

qxkcd <- function(p, sd = 1, log.p = FALSE, swap.end.points = FALSE){
  if(length(sd) == 1){
    sd = rep(sd, length(p))
  }
  if(length(sd) < length(p)){
    sd = append(sd, rep(1, length(length(p) - length(sd))))
  }
  
  out <- .Call("qxkcd", r_p = as.double(p), r_d = as.double(sd), r_log= as.integer(log.p), r_swap_end_points= as.integer(swap.end.points)
      , r_n = length(p), r_sd_n = length(sd), PACKAGE = "xkcd")
  return(out)
}



