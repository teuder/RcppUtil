df <- data.frame(matrix(rep(LETTERS,10000),26))

nrow <- 1000000
ncol <- 100
y <- sample(c(0,1),  nrow, replace = TRUE)
df <- data.frame(y, matrix(sample(1:5,  nrow*ncol, replace = TRUE), nrow=nrow, ncol=ncol))
res2 <- apply_catdap(y,df)

res <- catdap::catdap1(df,"y",FALSE)$aic[1,]

microbenchmark::microbenchmark(catdap::catdap1(df,"y",FALSE)$aic[1,],apply_catdap(y,df))
