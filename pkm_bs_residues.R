# Residue numbers are converted to use a zero-based index.
# The resShift value is subtracted so the first residue is 0.
pyr = function(r, bSort = FALSE, bChain = FALSE)
{
  resShift = 13
  res = list()
  res[[1]]  = c()
  res[[2]]  = c()
  res[[3]]  = c(269, 271, 294, 295, 327) - resShift
  res[[4]]  = c(290, 292, 359) - resShift
  res[[5]]  = c(72, 361) - resShift
  res[[6]]  = c(112, 242) - resShift
  res[[7]]  = c(74, 117, 243, 328) - resShift
  res[[8]]  = c(76, 114, 209, 298) - resShift
  res[[9]]  = c(177, 362) - resShift
  res[[10]] = c(119, 363) - resShift
  res[[11]] = c(77, 207, 331, 365) - resShift
  res[[12]] = c(176, 366) - resShift
  res[[13]] = c(206) - resShift
  
  if (bChain)
  {
    res[[7]] = c(res[[7]], 341 - resShift)
  }
  
  i   = 1
  ret = c()
  while (i <= r)
  {
    ret = c(ret, res[[i]])
    i   = i + 1
  }
  
  if (bSort) {ret = sort(ret)}
  
  return(ret)
}


fbp = function(r, bSort = FALSE)
{
  resShift = 14
  res = list()
  res[[1]]  = c()
  res[[2]]  = c()
  res[[3]]  = c(432, 437, 482, 514, 516, 518, 519, 520, 521) - resShift
  res[[4]]  = c(431, 433, 434, 436, 489, 517, 522) - resShift
  res[[5]]  = c(513) - resShift
  res[[6]]  = c(409, 455, 486) - resShift
  res[[7]]  = c(440, 453) - resShift
  res[[8]]  = c(454) - resShift
  res[[9]]  = c() - resShift
  res[[10]] = c() - resShift
  res[[11]] = c() - resShift
  res[[12]] = c() - resShift
  res[[13]] = c() - resShift
  
  i   = 1
  ret = c()
  while (i <= r)
  {
    ret = c(ret, res[[i]])
    i   = i + 1
  }
  
  if (bSort) {ret = sort(ret)}
  
  return(ret)
}


oxl = function(r, bSort = FALSE, bChain = FALSE)
{
  resShift = 14
  res = list()
  res[[1]]  = c()
  res[[2]]  = c()
  res[[3]]  = c(270, 272, 296, 328) - resShift
  res[[4]]  = c(293, 295) - resShift
  res[[5]]  = c(73, 291, 327, 360) - resShift
  res[[6]]  = c(362) - resShift
  res[[7]]  = c(114, 178) - resShift
  res[[8]]  = c(75, 77, 115, 118) - resShift
  res[[9]]  = c() - resShift
  res[[10]] = c() - resShift
  res[[11]] = c() - resShift
  res[[12]] = c() - resShift
  res[[13]] = c() - resShift
  
  if (bChain)
  {
    res[[7]] = c(res[[7]], 342 - resShift)
  }
  
  i   = 1
  ret = c()
  while (i <= r)
  {
    ret = c(ret, res[[i]])
    i   = i + 1
  }
  
  if (bSort) {ret = sort(ret)}
  
  return(ret)
}


gol = function(r, bSort = FALSE, bChain = FALSE)
{
  resShift = 14
  res = list()
  res[[1]]  = c()
  res[[2]]  = c()
  res[[3]]  = c(469) - resShift
  res[[4]]  = c(43, 44, 46, 70, 464, 466, 470, 471) - resShift
  res[[5]]  = c(468) - resShift
  res[[6]]  = c(68) - resShift
  res[[7]]  = c(106, 107, 379) - resShift
  res[[8]]  = c(500) - resShift
  res[[9]]  = c() - resShift
  res[[10]] = c() - resShift
  res[[11]] = c() - resShift
  res[[12]] = c() - resShift
  res[[13]] = c() - resShift
  
  i   = 1
  ret = c()
  while (i <= r)
  {
    ret = c(ret, res[[i]])
    i   = i + 1
  }
  
  if (bSort) {ret = sort(ret)}
  
  return(ret)
}


ser = function(r, bSort = FALSE, bChain = FALSE)
{
  resShift = 14
  res = list()
  res[[1]]  = c()
  res[[2]]  = c()
  res[[3]]  = c(43, 44, 70, 106, 464, 469) - resShift
  res[[4]]  = c(46, 468, 470) - resShift
  res[[5]]  = c(466) - resShift
  res[[6]]  = c(68) - resShift
  res[[7]]  = c() - resShift
  res[[8]]  = c() - resShift
  res[[9]]  = c() - resShift
  res[[10]] = c() - resShift
  res[[11]] = c() - resShift
  res[[12]] = c() - resShift
  res[[13]] = c() - resShift
  
  i   = 1
  ret = c()
  while (i <= r)
  {
    ret = c(ret, res[[i]])
    i   = i + 1
  }
  
  if (bSort) {ret = sort(ret)}
  
  return(ret)
}


po4 = function(r, bSort = FALSE, bChain = FALSE)
{
  resShift = 13
  res = list()
  res[[1]]  = c()
  res[[2]]  = c()
  res[[3]]  = c(431, 432, 433, 435, 436) - resShift
  res[[4]]  = c(519) - resShift
  res[[5]]  = c(430, 520, 521) - resShift
  res[[6]]  = c() - resShift
  res[[7]]  = c(513) - resShift
  res[[8]]  = c(481, 516, 518) - resShift
  res[[9]]  = c() - resShift
  res[[10]] = c() - resShift
  res[[11]] = c() - resShift
  res[[12]] = c() - resShift
  res[[13]] = c() - resShift
  
  i   = 1
  ret = c()
  while (i <= r)
  {
    ret = c(ret, res[[i]])
    i   = i + 1
  }
  
  if (bSort) {ret = sort(ret)}
  
  return(ret)
}


atp = function(r, bSort = FALSE, bChain = FALSE)
{
  resShift = 14
  res = list()
  res[[1]]  = c()
  res[[2]]  = c()
  res[[3]]  = c(75, 78, 120, 207) - resShift
  res[[4]]  = c(50, 51, 53, 73, 83, 84, 178, 270, 296, 362, 363, 366, 367) - resShift
  res[[5]]  = c(77, 79, 113, 128, 177, 205, 243, 272) - resShift
  res[[6]]  = c(114, 118, 129, 204, 208, 244, 295, 360) - resShift
  res[[7]]  = c(131, 293, 335, 364) - resShift
  res[[8]]  = c(123, 130, 291, 294, 327, 329, 332) - resShift
  res[[9]]  = c() - resShift
  res[[10]] = c() - resShift
  res[[11]] = c() - resShift
  res[[12]] = c() - resShift
  res[[13]] = c() - resShift
  
  if (bChain)
  {
    res[[10]] = c(res[[10]], 342 - resShift)
  }
  
  i   = 1
  ret = c()
  while (i <= r)
  {
    ret = c(ret, res[[i]])
    i   = i + 1
  }
  
  if (bSort) {ret = sort(ret)}
  
  return(ret)
}
