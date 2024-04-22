popnormINDEX <- function (age = NA, sex = "B", region = "BE", year = 2018) 
{
  
  #be_qol <- popnormINDEX(age = 15:100,sex='B',region='BE',year=2018)
  #saveRDS(be_qol,'data/eq5d/popnormINDEX.rds')
  be_qol <- readRDS('data/eq5d/popnormINDEX.rds')
  return(be_qol)
}

